program simple_test_continuous_3D_pcg_reconstruction
use continuous_3D_pcg_reconstruction_volume_test, only: run_volume_fixture
use continuous_3D_pcg_reconstruction_noise_test, only: run_volume_noise
use continuous_3D_pcg_reconstruction_halfset_test, only: run_halfset_fsc
use continuous_3D_pcg_reconstruction_test_helpers, only: run_current_executable_case
implicit none
character(len=256) :: argument, selected_case
integer :: case_count, failures, groups_passed, groups_run, groups_skipped, iarg, separator

selected_case = ''
case_count = 0
do iarg = 1, command_argument_count()
    call get_command_argument(iarg,argument)
    separator = index(argument,'=')
    if( separator > 0 )then
        if( trim(argument(:separator-1)) == 'case' )then
            case_count = case_count+1
            selected_case = trim(argument(separator+1:))
        endif
    endif
enddo
if( case_count > 1 ) &
    &error stop 'continuous 3D PCG reconstruction suite accepts only one case= argument'
if( case_count == 1 .and. len_trim(selected_case) == 0 ) &
    &error stop 'continuous 3D PCG reconstruction suite requires a nonempty case= value'
if( case_count == 0 )then
    failures = 0
    groups_passed = 0
    groups_run = 0
    groups_skipped = 0
    write(*,'(a)') 'Continuous 3D PCG reconstruction test suite: START'
    call run_current_executable_case('volume_fixture',groups_run,groups_passed,groups_skipped,failures)
    call run_current_executable_case('volume_noise',groups_run,groups_passed,groups_skipped,failures)
    call run_current_executable_case('halfset_fsc',groups_run,groups_passed,groups_skipped,failures)
    write(*,'(/,a)') 'Continuous 3D PCG reconstruction test suite: SUMMARY'
    write(*,'(a,i0)') '  Test groups scheduled: ',3
    write(*,'(a,i0)') '  Test groups run:       ',groups_run
    write(*,'(a,i0)') '  Test groups passed:    ',groups_passed
    write(*,'(a,i0)') '  Test groups skipped:   ',groups_skipped
    write(*,'(a,i0)') '  Test groups failed:    ',failures
    if( failures /= 0 ) error stop 1
    write(*,'(a)') 'Continuous 3D PCG reconstruction test suite: PASS'
else
    select case(trim(selected_case))
    case('volume_fixture')
        call run_volume_fixture()
    case('volume_noise')
        call run_volume_noise()
    case('halfset_fsc')
        call run_halfset_fsc()
    case default
        error stop 'unknown continuous 3D PCG reconstruction case: '//trim(selected_case)
    end select
endif
end program simple_test_continuous_3D_pcg_reconstruction
