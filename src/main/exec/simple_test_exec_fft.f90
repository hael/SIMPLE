!@descr: execution of test fft processing commanders
module simple_test_exec_fft
use simple_cmdline,             only: cmdline
use simple_commanders_test_fft, only: commander_test_corrs2weights_test, &
                                      commander_test_eval_polarftcc, commander_test_ft_expanded, &
                                      commander_test_gencorrs_fft, commander_test_order_corr, &
                                      commander_test_phasecorr, commander_test_polarops, &
                                      commander_test_rank_weights, commander_test_rotate_ref
implicit none

public :: exec_fft_commander
private

type(commander_test_corrs2weights_test) :: xcorrs2weights_test
type(commander_test_eval_polarftcc)     :: xeval_polarftcc
type(commander_test_ft_expanded)        :: xft_expanded
type(commander_test_gencorrs_fft)       :: xgencorrs_fft
type(commander_test_order_corr)         :: xorder_corr
type(commander_test_phasecorr)          :: xphasecorr
type(commander_test_polarops)           :: xpolarops
type(commander_test_rank_weights)       :: xrank_weights
type(commander_test_rotate_ref)         :: xrotate_ref

contains

    subroutine exec_fft_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'corrs2weights_test' )
                call xcorrs2weights_test%execute(cline)
            case( 'eval_polarftcc' )
                call xeval_polarftcc%execute(cline)
            case( 'ft_expanded' )
                call xft_expanded%execute(cline)
            case( 'gencorrs_fft' )
                call xgencorrs_fft%execute(cline)
            case( 'order_corr' )
                call xorder_corr%execute(cline)
            case( 'phasecorr' )
                call xphasecorr%execute(cline)
            case( 'polarops' )
                call xpolarops%execute(cline)
            case( 'rank_weights' )
                call xrank_weights%execute(cline)
            case( 'rotate_ref' )
                call xrotate_ref%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_fft_commander

end module simple_test_exec_fft
