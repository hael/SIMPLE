!@descr: module defining the user interfaces for fft testprograms in the simple_test_exec suite
module simple_test_ui_fft
use simple_ui_modules
implicit none

type(ui_program), target :: corrs2weights_test
type(ui_program), target :: eval_polarftcc
type(ui_program), target :: ft_expanded
type(ui_program), target :: gencorrs_fft
type(ui_program), target :: order_corr
type(ui_program), target :: phasecorr
type(ui_program), target :: polarops
type(ui_program), target :: rank_weights
type(ui_program), target :: rotate_ref

contains

    subroutine construct_test_fft_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_corrs2weights_test(tsttab)
        call new_eval_polarftcc(tsttab)
        call new_ft_expanded(tsttab)
        call new_gencorrs_fft(tsttab)
        call new_order_corr(tsttab)
        call new_phasecorr(tsttab)
        call new_polarops(tsttab)
        call new_rank_weights(tsttab)
        call new_rotate_ref(tsttab)
    end subroutine construct_test_fft_programs

    subroutine print_test_fft_programs( logfhandle )
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('FFT:', C_UNDERLINED)
        write(logfhandle,'(A)') corrs2weights_test%name%to_char()
        write(logfhandle,'(A)') eval_polarftcc%name%to_char()
        write(logfhandle,'(A)') ft_expanded%name%to_char()
        write(logfhandle,'(A)') gencorrs_fft%name%to_char()
        write(logfhandle,'(A)') order_corr%name%to_char()
        write(logfhandle,'(A)') phasecorr%name%to_char()
        write(logfhandle,'(A)') polarops%name%to_char()
        write(logfhandle,'(A)') rank_weights%name%to_char()
        write(logfhandle,'(A)') rotate_ref%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_fft_programs

    subroutine new_corrs2weights_test( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call corrs2weights_test%new(&
        &'corrs2weights_test',&                ! name
        &'corrs2weights_test ',&               ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call corrs2weights_test%add_input(UI_IO, )
        ! parameter input/output
        !call corrs2weights_test%add_input(UI_IMG, )
        ! alternative inputs
        !call corrs2weights_test%add_input(UI_PARM, )
        ! search controls
        !call corrs2weights_test%add_input(UI_SRCH, )
        ! filter controls
        !call corrs2weights_test%add_input(UI_FILT, )
        ! mask controls
        !call corrs2weights_test%add_input(UI_MASK, )
        ! computer controls
        !call corrs2weights_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('corrs2weights_test', corrs2weights_test, tsttab)
    end subroutine new_corrs2weights_test

    subroutine new_eval_polarftcc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call eval_polarftcc%new(&
        &'eval_polarftcc',&                    ! name
        &'eval_polarftcc ',&                   ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call eval_polarftcc%add_input(UI_IO, )
        ! parameter input/output
        !call eval_polarftcc%add_input(UI_IMG, )
        ! alternative inputs
        !call eval_polarftcc%add_input(UI_PARM, )
        ! search controls
        !call eval_polarftcc%add_input(UI_SRCH, )
        ! filter controls
        !call eval_polarftcc%add_input(UI_FILT, )
        ! mask controls
        !call eval_polarftcc%add_input(UI_MASK, )
        ! computer controls
        !call eval_polarftcc%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('eval_polarftcc', eval_polarftcc, tsttab)
    end subroutine new_eval_polarftcc

    subroutine new_ft_expanded( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call ft_expanded%new(&
        &'ft_expanded',&                       ! name
        &'ft_expanded ',&                      ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call ft_expanded%add_input(UI_IO, )
        ! parameter input/output
        !call ft_expanded%add_input(UI_IMG, )
        ! alternative inputs
        !call ft_expanded%add_input(UI_PARM, )
        ! search controls
        !call ft_expanded%add_input(UI_SRCH, )
        ! filter controls
        !call ft_expanded%add_input(UI_FILT, )
        ! mask controls
        !call ft_expanded%add_input(UI_MASK, )
        ! computer controls
        !call ft_expanded%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('ft_expanded', ft_expanded, tsttab)
    end subroutine new_ft_expanded

    subroutine new_gencorrs_fft( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call gencorrs_fft%new(&
        &'gencorrs_fft',&                      ! name
        &'gencorrs_fft ',&                     ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call gencorrs_fft%add_input(UI_IO, )
        ! parameter input/output
        !call gencorrs_fft%add_input(UI_IMG, )
        ! alternative inputs
        !call gencorrs_fft%add_input(UI_PARM, )
        ! search controls
        !call gencorrs_fft%add_input(UI_SRCH, )
        ! filter controls
        !call gencorrs_fft%add_input(UI_FILT, )
        ! mask controls
        !call gencorrs_fft%add_input(UI_MASK, )
        ! computer controls
        !call gencorrs_fft%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('gencorrs_fft', gencorrs_fft, tsttab)
    end subroutine new_gencorrs_fft

    subroutine new_order_corr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call order_corr%new(&
        &'order_corr',&                        ! name
        &'order_corr ',&                       ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call order_corr%add_input(UI_IO, )
        ! parameter input/output
        !call order_corr%add_input(UI_IMG, )
        ! alternative inputs
        !call order_corr%add_input(UI_PARM, )
        ! search controls
        !call order_corr%add_input(UI_SRCH, )
        ! filter controls
        !call order_corr%add_input(UI_FILT, )
        ! mask controls
        !call order_corr%add_input(UI_MASK, )
        ! computer controls
        !call order_corr%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('order_corr', order_corr, tsttab)
    end subroutine new_order_corr

    subroutine new_phasecorr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call phasecorr%new(&
        &'phasecorr',&                         ! name
        &'phasecorr ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call phasecorr%add_input(UI_IO, )
        ! parameter input/output
        !call phasecorr%add_input(UI_IMG, )
        ! alternative inputs
        !call phasecorr%add_input(UI_PARM, )
        ! search controls
        !call phasecorr%add_input(UI_SRCH, )
        ! filter controls
        !call phasecorr%add_input(UI_FILT, )
        ! mask controls
        !call phasecorr%add_input(UI_MASK, )
        ! computer controls
        !call phasecorr%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('phasecorr', phasecorr, tsttab)
    end subroutine new_phasecorr

    subroutine new_polarops( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call polarops%new(&
        &'polarops',&                         ! name
        &'polarops ',&                        ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                 ! executable
        &.false.)                             ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call polarops%add_input(UI_IO, )
        ! parameter input/output
        !call polarops%add_input(UI_IMG, )
        ! alternative inputs
        !call polarops%add_input(UI_PARM, )
        ! search controls
        !call polarops%add_input(UI_SRCH, )
        ! filter controls
        !call polarops%add_input(UI_FILT, )
        ! mask controls
        !call polarops%add_input(UI_MASK, )
        ! computer controls
        !call polarops%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('polarops', polarops, tsttab)
    end subroutine new_polarops

    subroutine new_rank_weights( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call rank_weights%new(&
        &'rank_weights',&                      ! name
        &'rank_weights ',&                     ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call rank_weights%add_input(UI_IO, )
        ! parameter input/output
        !call rank_weights%add_input(UI_IMG, )
        ! alternative inputs
        !call rank_weights%add_input(UI_PARM, )
        ! search controls
        !call rank_weights%add_input(UI_SRCH, )
        ! filter controls
        !call rank_weights%add_input(UI_FILT, )
        ! mask controls
        !call rank_weights%add_input(UI_MASK, )
        ! computer controls
        !call rank_weights%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('rank_weights', rank_weights, tsttab)
    end subroutine new_rank_weights

    subroutine new_rotate_ref( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call rotate_ref%new(&
        &'rotate_ref',&                        ! name
        &'rotate_ref ',&                       ! descr_short
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call rotate_ref%add_input(UI_IO, )
        ! parameter input/output
        !call rotate_ref%add_input(UI_IMG, )
        ! alternative inputs
        !call rotate_ref%add_input(UI_PARM, )
        ! search controls
        !call rotate_ref%add_input(UI_SRCH, )
        ! filter controls
        !call rotate_ref%add_input(UI_FILT, )
        ! mask controls
        !call rotate_ref%add_input(UI_MASK, )
        ! computer controls
        !call rotate_ref%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('rotate_ref', rotate_ref, tsttab)
    end subroutine new_rotate_ref

end module simple_test_ui_fft
