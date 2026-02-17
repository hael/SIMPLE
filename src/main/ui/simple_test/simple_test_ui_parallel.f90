!@ descr: module defining the user interfaces for parallel programs in the simple_test_exec suite
module simple_test_ui_parallel
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_openmp
type(ui_program), target :: simple_test_openacc
type(ui_program), target :: simple_test_coarrays
type(ui_program), target :: simple_test_simd

contains

    subroutine construct_parallel_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_openmp(prgtab)
        call new_simple_test_openacc(prgtab)
        call new_simple_test_coarrays(prgtab)
        call new_simple_test_simd(prgtab)
    end subroutine construct_parallel_programs

    subroutine print_parallel_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('PARALLEL:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_openmp%name%to_char()
        write(logfhandle,'(A)') simple_test_openacc%name%to_char()
        write(logfhandle,'(A)') simple_test_coarrays%name%to_char()
        write(logfhandle,'(A)') simple_test_simd%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_parallel_programs

    subroutine new_simple_test_openmp( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_openmp', simple_test_openmp, prgtab)
    end subroutine new_simple_test_openmp

    subroutine new_simple_test_openacc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_openacc', simple_test_openacc, prgtab)
    end subroutine new_simple_test_openacc

    subroutine new_simple_test_coarrays( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_coarrays', simple_test_coarrays, prgtab)
    end subroutine new_simple_test_coarrays

    subroutine new_simple_test_simd( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_simd', simple_test_simd, prgtab)
    end subroutine new_simple_test_simd

end module simple_test_ui_parallel
