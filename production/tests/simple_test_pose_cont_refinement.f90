program simple_test_pose_cont_refinement
    use pose_cont_refinement_numerics_test, only: run_pose_cont_numerics
    use pose_cont_refinement_solver_test,   only: run_pose_cont_solver
    implicit none

    character(len=32) :: selected_case
    integer :: occurrences

    call find_selected_case(selected_case,occurrences)
    if( occurrences > 1 ) error stop 'pose-cont suite accepts only one case= argument'
    if( occurrences == 0 )then
        call run_pose_cont_numerics()
        call run_pose_cont_solver()
        write(*,'(a)') 'POSE_CONT_REFINEMENT_SUITE: PASS'
    else
        select case(trim(selected_case))
        case('numerics')
            call run_pose_cont_numerics()
        case('solver')
            call run_pose_cont_solver()
        case default
            error stop 'pose-cont suite requires case=numerics or case=solver'
        end select
    endif

contains

    subroutine find_selected_case(case_name,count)
        character(len=*), intent(out) :: case_name
        integer, intent(out) :: count
        character(len=256) :: argument
        integer :: iarg, separator, status

        case_name = ''
        count = 0
        do iarg = 1, command_argument_count()
            call get_command_argument(iarg,argument,status=status)
            if( status /= 0 ) error stop 'could not read pose-cont test argument'
            separator = index(argument,'=')
            if( separator <= 1 ) cycle
            if( trim(argument(:separator-1)) /= 'case' ) cycle
            count = count+1
            case_name = trim(argument(separator+1:))
        enddo
    end subroutine find_selected_case

end program simple_test_pose_cont_refinement
