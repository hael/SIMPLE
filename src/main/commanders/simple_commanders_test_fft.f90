!@descr: for all fft tests
module simple_commanders_test_fft
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_corrs2weights_test
  contains
    procedure :: execute      => exec_test_corrs2weights_test
end type commander_test_corrs2weights_test

type, extends(commander_base) :: commander_test_eval_polarftcc
  contains
    procedure :: execute      => exec_test_eval_polarftcc
end type commander_test_eval_polarftcc

type, extends(commander_base) :: commander_test_ft_expanded
  contains
    procedure :: execute      => exec_test_ft_expanded
end type commander_test_ft_expanded

type, extends(commander_base) :: commander_test_gencorrs_fft
  contains
    procedure :: execute      => exec_test_gencorrs_fft
end type commander_test_gencorrs_fft

type, extends(commander_base) :: commander_test_order_corr
  contains
    procedure :: execute      => exec_test_order_corr
end type commander_test_order_corr

type, extends(commander_base) :: commander_test_phasecorr
  contains
    procedure :: execute      => exec_test_phasecorr
end type commander_test_phasecorr

type, extends(commander_base) :: commander_test_polarops
  contains
    procedure :: execute      => exec_test_polarops
end type commander_test_polarops

type, extends(commander_base) :: commander_test_rank_weights
  contains
    procedure :: execute      => exec_test_rank_weights
end type commander_test_rank_weights

type, extends(commander_base) :: commander_test_rotate_ref
  contains
    procedure :: execute      => exec_test_rotate_ref
end type commander_test_rotate_ref

contains

    subroutine exec_test_corrs2weights_test( self, cline )
        use simple_core_module_api
        class(commander_test_corrs2weights_test), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        real    :: corrs(12), weights(12)
        integer :: i
        corrs(1)  = -1.
        corrs(2)  = 0.0
        corrs(3)  = 0.005
        corrs(4)  = 0.1
        corrs(5)  = 0.2
        corrs(6)  = 0.3
        corrs(7)  = 0.4
        corrs(8)  = 0.5
        corrs(9)  = 0.51
        corrs(10) = 0.52
        corrs(11) = 0.53
        corrs(12) = 0.6
        weights = corrs2weights(corrs, CORRW_CRIT)
        do i=1,size(corrs)
            print *, 'corr/weight: ', corrs(i), weights(i)
        end do
        call simple_end('**** SIMPLE_TEST_CORRS2WEIGHTS_TEST_WORKFLOW NORMAL STOP ****')
    end subroutine exec_test_corrs2weights_test

subroutine exec_test_eval_polarftcc( self, cline )
    class(commander_test_eval_polarftcc), intent(inout) :: self
    class(cmdline),                       intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_EVAL_POLARFTCC_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_eval_polarftcc

subroutine exec_test_ft_expanded( self, cline )
    class(commander_test_ft_expanded), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_FT_EXPANDED_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ft_expanded

subroutine exec_test_gencorrs_fft( self, cline )
    class(commander_test_gencorrs_fft), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_GENCORRS_FFT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_gencorrs_fft

subroutine exec_test_order_corr( self, cline )
    class(commander_test_order_corr), intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ORDER_CORR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_order_corr

subroutine exec_test_phasecorr( self, cline )
    class(commander_test_phasecorr), intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_PHASECORR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_phasecorr

subroutine exec_test_polarops( self, cline )
    class(commander_test_polarops), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_POLAROPS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_polarops

subroutine exec_test_rank_weights( self, cline )
    class(commander_test_rank_weights), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_RANK_WEIGHTS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_rank_weights

subroutine exec_test_rotate_ref( self, cline )
    class(commander_test_rotate_ref), intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ROTATE_REF_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_rotate_ref

end module simple_commanders_test_fft
