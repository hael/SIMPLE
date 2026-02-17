!@ descr: module defining the user interfaces for FFTprograms in the simple_test_exec suite
module simple_test_ui_fft
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_phasecorr
type(ui_program), target :: simple_test_order_corr
type(ui_program), target :: simple_test_gencorrs_fft
type(ui_program), target :: simple_test_ft_expanded
type(ui_program), target :: simple_test_eval_polarftcc
type(ui_program), target :: simple_test_polarops
type(ui_program), target :: simple_test_corrs2weights
type(ui_program), target :: simple_test_rank_weights
type(ui_program), target :: simple_test_rotate_ref

contains

    subroutine construct_fft_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_phasecorr(prgtab)
        call new_simple_test_order_corr(prgtab)
        call new_simple_test_gencorrs_fft(prgtab)
        call new_simple_test_ft_expanded(prgtab)
        call new_simple_test_eval_polarftcc(prgtab)
        call new_simple_test_polarops(prgtab)
        call new_simple_test_corrs2weights(prgtab)
        call new_simple_test_rank_weights(prgtab)
        call new_simple_test_rotate_ref(prgtab)
    end subroutine construct_fft_programs

    subroutine print_fft_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('FFT:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_phasecorr%name%to_char()
        write(logfhandle,'(A)') simple_test_order_corr%name%to_char()
        write(logfhandle,'(A)') simple_test_gencorrs_fft%name%to_char()
        write(logfhandle,'(A)') simple_test_ft_expanded%name%to_char()
        write(logfhandle,'(A)') simple_test_eval_polarftcc%name%to_char()
        write(logfhandle,'(A)') simple_test_polarops%name%to_char()
        write(logfhandle,'(A)') simple_test_corrs2weights%name%to_char()
        write(logfhandle,'(A)') simple_test_rank_weights%name%to_char()
        write(logfhandle,'(A)') simple_test_rotate_ref%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_fft_programs

    subroutine new_simple_test_phasecorr( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_phasecorr', simple_test_phasecorr, prgtab)
    end subroutine new_simple_test_phasecorr

    subroutine new_simple_test_order_corr( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_order_corr', simple_test_order_corr, prgtab)
    end subroutine new_simple_test_order_corr

    subroutine new_simple_test_gencorrs_fft( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_gencorrs_fft', simple_test_gencorrs_fft, prgtab)
    end subroutine new_simple_test_gencorrs_fft

    subroutine new_simple_test_ft_expanded( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_ft_expanded', simple_test_ft_expanded, prgtab)
    end subroutine new_simple_test_ft_expanded

    subroutine new_simple_test_eval_polarftcc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_eval_polarftcc', simple_test_eval_polarftcc, prgtab)
    end subroutine new_simple_test_eval_polarftcc

    subroutine new_simple_test_polarops( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_polarops', simple_test_polarops, prgtab)
    end subroutine new_simple_test_polarops

    subroutine new_simple_test_corrs2weights( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_corrs2weights', simple_test_corrs2weights, prgtab)
    end subroutine new_simple_test_corrs2weights

    subroutine new_simple_test_rank_weights( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_rank_weights', simple_test_rank_weights, prgtab)
    end subroutine new_simple_test_rank_weights

    subroutine new_simple_test_rotate_ref( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_rotate_ref', simple_test_rotate_ref, prgtab)
    end subroutine new_simple_test_rotate_ref

end module simple_test_ui_fft
