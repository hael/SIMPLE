!@descr: execution of test parallel processing commanders
module simple_test_exec_parallel
use simple_cmdline,                  only: cmdline
use simple_commanders_test_parallel, only: commander_test_coarrays, commander_test_openacc, &
                                           commander_test_openmp, commander_test_simd
implicit none

public :: exec_parallel_commander
private

type(commander_test_coarrays) :: xcoarrays
type(commander_test_openacc)  :: xopenacc
type(commander_test_openmp)   :: xopenmp
type(commander_test_simd)     :: xsimd

contains

    subroutine exec_parallel_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'coarrays' )
                call xcoarrays%execute(cline)
            case( 'openacc' )
                call xopenacc%execute(cline)
            case( 'openmp' )
                call xopenmp%execute(cline)
            case( 'simd' )
                call xsimd%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_parallel_commander

end module simple_test_exec_parallel
