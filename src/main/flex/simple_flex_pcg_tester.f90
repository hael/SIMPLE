!@descr: unit and library tests of the flex_pca PCG M-step operator (simple_flex_pca_pcg)
! Drives the module's white-box self-test test_flex_pcg_operator (it needs the private
! components of flex_pcg_t and the private scatter/fold kernels, so it stays in the
! module) and asserts its four checks by name: (A) the operator through the doubled-
! coordinate kernels against the exact nonuniform-DFT Gram, (B) the right-hand-side
! deposit against the exact adjoint, (C) at box <= 32 a preconditioned CG solve that
! recovers the volume, (D) the band-list kernels against the dense fold. The fast gate
! runs box 32 with 200 samples and the clean baseline solve; nightly runs box 64 with
! 400 samples (no solve) and the twelve-setting solve sweep at box 32, where every clean
! solve is held to the baseline criterion.
module simple_flex_pcg_tester
use simple_flex_pca_pcg, only: test_flex_pcg_operator
use simple_test_utils
implicit none
private
public :: run_all_flex_pcg_tests, run_all_flex_pcg_lib_tests, run_all_flex_pcg_sweep_tests

contains

    subroutine run_all_flex_pcg_tests()
        write(*,'(A)') '**** running all flex PCG operator tests ****'
        write(*,'(A)') 'test_flex_pcg_operator_box32'
        call check_operator(32, 200, .true., .false.)
    end subroutine run_all_flex_pcg_tests

    subroutine run_all_flex_pcg_lib_tests()
        write(*,'(A)') '**** running all flex PCG operator library tests ****'
        write(*,'(A)') 'test_flex_pcg_operator_box64'
        call check_operator(64, 400, .false., .false.)
    end subroutine run_all_flex_pcg_lib_tests

    subroutine run_all_flex_pcg_sweep_tests()
        write(*,'(A)') '**** running all flex PCG solve sweep tests ****'
        write(*,'(A)') 'test_flex_pcg_operator_box32_sweep'
        call check_operator(32, 200, .true., .true.)
    end subroutine run_all_flex_pcg_sweep_tests

    subroutine check_operator( box, nsamples, with_solve, sweep )
        integer, intent(in) :: box, nsamples
        logical, intent(in) :: with_solve, sweep
        logical :: l_pass, passes(4)
        character(len=32) :: tag
        write(tag,'(A,I0,A,I0,A)') 'box ', box, ', ', nsamples, ' samples: '
        call test_flex_pcg_operator(box, nsamples, l_pass, passes=passes, sweep=sweep)
        call assert_true(passes(1), trim(tag)//' (A) kernel operator matches the exact Gram (scale 5%, residual 10%)')
        call assert_true(passes(2), trim(tag)//' (B) rhs deposit matches the exact adjoint (scale 5%, residual 10%)')
        if( with_solve )then
            if( sweep )then
                call assert_true(passes(3), trim(tag)//' (C) every clean solve of the sweep recovers the volume on its support')
            else
                call assert_true(passes(3), trim(tag)//' (C) the CG solve recovers the volume on its support')
            endif
        endif
        call assert_true(passes(4), trim(tag)//' (D) band-list kernels and rhs equal the dense fold')
        call assert_true(l_pass .eqv. all(passes), trim(tag)//' the overall verdict is the conjunction of the four')
    end subroutine check_operator

end module simple_flex_pcg_tester
