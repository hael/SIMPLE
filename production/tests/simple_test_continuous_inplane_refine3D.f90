program simple_test_continuous_inplane_refine3D
use continuous_inplane_refine3D_baseline_test, only: run_legacy_baseline
use continuous_inplane_refine3D_policy_test, only: run_policy_gate, run_policy_rejection
use continuous_inplane_refine3D_state_test, only: run_search_state_contract
use continuous_inplane_refine3D_joint_test, only: run_joint_state_contract
use continuous_inplane_refine3D_direct_test, only: run_direct_route_contract
use continuous_inplane_refine3D_metadata_test, only: run_metadata_state_contract, &
    &run_metadata_project_contract
use continuous_inplane_refine3D_recovery_test, only: run_synthetic_recovery_contract
implicit none

character(len=4096) :: argument
character(len=256)  :: selected_case
integer :: iarg, separator

selected_case = ''
do iarg = 1, command_argument_count()
    call get_command_argument(iarg, argument)
    separator = index(argument, '=')
    if( separator > 1 )then
        if( trim(argument(:separator-1)) == 'case' ) &
            &selected_case = trim(argument(separator+1:))
    endif
enddo

if( len_trim(selected_case) > 0 )then
    call run_case(selected_case)
else
    call run_suite()
endif

contains

    subroutine run_suite()
        integer :: failures, groups_run, groups_skipped

        failures = 0
        groups_run = 0
        groups_skipped = 0
        write(*,'(a)') 'Continuous in-plane refine3D regression suite: START'
        if( argument_defined('projfile') .and. &
            &argument_defined('expected_snapshot') )then
            call run_case_in_child('baseline', failures)
            groups_run = groups_run + 1
        else
            groups_skipped = groups_skipped + 1
            write(*,'(/,a)') '>>> TEST [baseline] SKIP'
            write(*,'(a)') '    NOTICE: baseline comparison requires projfile and expected_snapshot.'
            write(*,'(a)') '    Generate the expected TSV first with case=baseline and snapshot=xx.tsv.'
            write(*,'(a)') '    The self-contained policy and search tests will continue.'
        endif
        call run_case_in_child('policy_suite', failures)
        groups_run = groups_run + 1
        call run_case_in_child('search_state', failures)
        groups_run = groups_run + 1
        call run_case_in_child('joint_state', failures)
        groups_run = groups_run + 1
        call run_case_in_child('direct_route', failures)
        groups_run = groups_run + 1
        call run_case_in_child('metadata_state', failures)
        groups_run = groups_run + 1
        if( argument_defined('vol1') )then
            call run_case_in_child('synthetic_recovery', failures)
            groups_run = groups_run + 1
        else
            groups_skipped = groups_skipped + 1
            write(*,'(/,a)') '>>> TEST [synthetic_recovery] SKIP'
            write(*,'(a)') &
                &'    NOTICE: synthetic recovery requires vol1=/absolute/path/reference_volume.mrc.'
            write(*,'(a)') '    The remaining self-contained tests will continue.'
        endif
        write(*,'(/,a)') 'Continuous in-plane refine3D regression suite: SUMMARY'
        write(*,'(a,i0)') '  Test groups scheduled: ', 7
        write(*,'(a,i0)') '  Test groups run:       ', groups_run
        write(*,'(a,i0)') '  Test groups skipped:   ', groups_skipped
        write(*,'(a,i0)') '  Test groups passed:    ', groups_run - failures
        write(*,'(a,i0)') '  Tests failed: ', failures
        if( failures /= 0 )then
            write(*,'(a)') 'Continuous in-plane refine3D regression suite: FAIL'
            error stop 1
        endif
        write(*,'(a)') 'Continuous in-plane refine3D regression suite: PASS'
    end subroutine run_suite

    logical function argument_defined(name) result(defined)
        character(len=*), intent(in) :: name
        character(len=4096) :: current_argument
        integer :: current_index, current_separator

        defined = .false.
        do current_index = 1, command_argument_count()
            call get_command_argument(current_index, current_argument)
            current_separator = index(current_argument, '=')
            if( current_separator <= 1 ) cycle
            if( trim(current_argument(:current_separator-1)) == trim(name) .and. &
                &len_trim(current_argument(current_separator+1:)) > 0 )then
                defined = .true.
                return
            endif
        enddo
    end function argument_defined

    subroutine run_case_in_child(label, failures)
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        character(len=4096) :: executable, child_command, parent_argument
        integer :: executable_len, command_status, exit_status, i

        call get_command_argument(0, executable, executable_len, command_status)
        if( command_status /= 0 .or. executable_len < 1 ) &
            &error stop 'could not determine refine3D suite executable path'
        child_command = '"'//trim(executable(:executable_len))//'" case='//trim(label)
        do i = 1, command_argument_count()
            call get_command_argument(i, parent_argument)
            if( index(parent_argument, 'case=') == 1 ) cycle
            child_command = trim(child_command)//' "'//trim(parent_argument)//'"'
        enddo
        write(*,'(/,a)') '>>> TEST ['//trim(label)//'] START'
        call execute_command_line(trim(child_command), wait=.true., &
            &exitstat=exit_status, cmdstat=command_status)
        if( command_status == 0 .and. exit_status == 0 )then
            write(*,'(a)') '>>> TEST ['//trim(label)//'] PASS'
        else
            failures = failures + 1
            write(*,'(a,i0,a,i0,a)') '>>> TEST ['//trim(label)// &
                &'] FAIL (cmdstat=', command_status, ', exitstat=', exit_status, ')'
        endif
    end subroutine run_case_in_child

    subroutine run_policy_suite()
        integer :: failures

        failures = 0
        write(*,'(a)') 'Continuous in-plane refine3D policy suite: START'
        call run_case_in_child('policy', failures)
        call run_rejection_in_child('policy_bad_value', 'inpl_cont must be yes or no', failures)
        call run_rejection_in_child('policy_hybrid', 'does not support objfun_den=yes', failures)
        call run_rejection_in_child('policy_denoised', 'requires ptcl_src=raw', failures)
        call run_rejection_in_child('policy_projrec', 'does not support projrec=yes', failures)
        if( failures /= 0 )then
            write(*,'(a,i0,a)') 'Continuous in-plane refine3D policy suite: FAIL (', &
                &failures, ' failure(s))'
            error stop 1
        endif
        write(*,'(a)') 'Continuous in-plane refine3D policy suite: PASS'
    end subroutine run_policy_suite

    subroutine run_rejection_in_child(label, expected_text, failures)
        use simple_syslib, only: get_process_id
        character(len=*), intent(in) :: label, expected_text
        integer, intent(inout) :: failures
        character(len=4096) :: executable, child_command, logfile
        integer :: executable_len, command_status, exit_status

        call get_command_argument(0, executable, executable_len, command_status)
        if( command_status /= 0 .or. executable_len < 1 ) &
            &error stop 'could not determine refine3D suite executable path'
        write(logfile,'(a,i0,a)') 'simple_test_continuous_inplane_refine3D_', &
            &get_process_id(), '_'//trim(label)//'.log'
        child_command = '"'//trim(executable(:executable_len))//'" case='//trim(label)// &
            &' > "'//trim(logfile)//'" 2>&1'
        write(*,'(/,a)') '>>> TEST ['//trim(label)//'] START'
        call execute_command_line(trim(child_command), wait=.true., &
            &exitstat=exit_status, cmdstat=command_status)
        if( command_status == 0 .and. exit_status /= 0 .and. &
            &file_contains(logfile, expected_text) )then
            write(*,'(a)') '>>> TEST ['//trim(label)//'] PASS'
        else
            failures = failures + 1
            write(*,'(a,i0,a,i0,a)') '>>> TEST ['//trim(label)// &
                &'] FAIL (cmdstat=', command_status, ', exitstat=', exit_status, ')'
        endif
    end subroutine run_rejection_in_child

    logical function file_contains(filename, expected_text) result(found)
        character(len=*), intent(in) :: filename, expected_text
        character(len=4096) :: line
        integer :: file_unit, ios

        found = .false.
        open(newunit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(file_unit,'(a)',iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, trim(expected_text)) > 0 )then
                found = .true.
                exit
            endif
        enddo
        close(file_unit)
    end function file_contains

    subroutine run_case(label)
        character(len=*), intent(in) :: label

        write(*,'(a)') '>>> Continuous in-plane refine3D child case: '//trim(label)
        select case(trim(label))
        case('baseline')
            if( .not. argument_defined('projfile') )then
                write(*,'(a)') 'REFINE3D_LEGACY_BASELINE SKIP: projfile was not supplied.'
                write(*,'(a)') 'NOTICE: rerun with projfile=/absolute/path/refine3D_output.simple.'
                return
            endif
            call run_legacy_baseline()
        case('policy')
            call run_policy_gate()
        case('policy_suite')
            call run_policy_suite()
        case('search_state')
            call run_search_state_contract()
        case('joint_state')
            call run_joint_state_contract()
        case('direct_route')
            call run_direct_route_contract()
        case('metadata_state')
            call run_metadata_state_contract()
        case('metadata_project')
            call run_metadata_project_contract()
        case('synthetic_recovery')
            call run_synthetic_recovery_contract()
        case('policy_bad_value', 'policy_hybrid', &
            &'policy_denoised', 'policy_projrec', 'policy_eval', &
            &'policy_sigma', 'policy_probabilistic', 'policy_neigh')
            call run_policy_rejection(trim(label))
        case default
            error stop 'unknown continuous in-plane refine3D regression case: '//trim(label)
        end select
    end subroutine run_case

end program simple_test_continuous_inplane_refine3D
