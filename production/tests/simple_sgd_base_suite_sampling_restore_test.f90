module base_sgd_sampling_restore_test
contains
    subroutine run_sampling_and_restore_policy()
        use simple_classaverager, only: cavger_zero_support_recovery
        use simple_oris, only: oris
        implicit none
        type(oris) :: os
        integer, allocatable :: inds(:)
        integer :: nsamp
        real :: frac

        call os%new(10, .true.)
        call os%set_all2single('state', 1.0)
        frac = 0.6
        call os%sample4update_rnd([1, 10], frac, nsamp, inds, .true.)
        if( nsamp /= 6 ) error stop 'SGD mini-batch fraction regression failed'
        if( any(inds < 1) .or. any(inds > 10) ) error stop 'SGD mini-batch index regression failed'
        call os%kill()

        if( .not. cavger_zero_support_recovery(0, .true.) ) &
            error stop 'zero-support recovery regression failed'
        if( cavger_zero_support_recovery(0, .false.) ) &
            error stop 'zero-support first-restore regression failed'
        if( cavger_zero_support_recovery(1, .true.) ) &
            error stop 'supported-class restoration regression failed'
        write(*,'(A)') 'SGD mini-batch and zero-support policy: PASS'
    end subroutine run_sampling_and_restore_policy
end module base_sgd_sampling_restore_test
