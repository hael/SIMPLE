!@descr: module defining the user interfaces for fft testprograms in the simple_test_exec suite
module simple_test_ui_fft
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('fft', 'FFT', 20)
type(ui_program), target :: corrs2weights_test
type(ui_program), target :: eval_polarftcc
type(ui_program), target :: ft_expanded
type(ui_program), target :: gencorrs_fft
type(ui_program), target :: order_corr
type(ui_program), target :: phasecorr
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
        call new_rank_weights(tsttab)
        call new_rotate_ref(tsttab)
    end subroutine construct_test_fft_programs

subroutine new_corrs2weights_test( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call corrs2weights_test%new(&
        &'corrs2weights_test',&                ! name
        &'corrs2weights test ',&               ! summary
        &'is a test program for generating correlation-based weights ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        call add_ui_program('corrs2weights_test', corrs2weights_test, tsttab, UI_CATEGORY)
    end subroutine new_corrs2weights_test

    subroutine new_eval_polarftcc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call eval_polarftcc%new(&
        &'eval_polarftcc',&                    ! name
        &'eval_polarftcc ',&                   ! summary
        &'is a test program for evaluating polar fourier cross-correlations ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call eval_polarftcc%add_input(UI_IO, )
        ! parameter input/output
        !call eval_polarftcc%add_input(UI_IMG, )
        ! <no additional inputs>
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
        call add_ui_program('eval_polarftcc', eval_polarftcc, tsttab, UI_CATEGORY)
    end subroutine new_eval_polarftcc

    subroutine new_ft_expanded( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call ft_expanded%new(&
        &'ft_expanded',&                       ! name
        &'ft_expanded ',&                      ! summary
        &'is a test program for shift search with L-BFGS-B using expanded Fourier transforms',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call ft_expanded%add_input(UI_IO, )
        ! parameter input/output
        !call ft_expanded%add_input(UI_IMG, )
        ! <no additional inputs>
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
        call add_ui_program('ft_expanded', ft_expanded, tsttab, UI_CATEGORY)
    end subroutine new_ft_expanded

    subroutine new_gencorrs_fft( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call gencorrs_fft%new(&
        &'gencorrs_fft',&                      ! name
        &'gencorrs_fft ',&                     ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call gencorrs_fft%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        call gencorrs_fft%add_input(UI_PARM, smpd, required_override=.true.)
        ! alternative inputs
        ! <no additional inputs>
        !call gencorrs_fft%add_input(UI_PARM, )
        ! search controls
        !call gencorrs_fft%add_input(UI_SRCH, )
        ! filter controls
        !call gencorrs_fft%add_input(UI_FILT, )
        ! mask controls
        call gencorrs_fft%add_input(UI_MASK, mskdiam, required_override=.true.)
        ! computer controls
        !call gencorrs_fft%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('gencorrs_fft', gencorrs_fft, tsttab, UI_CATEGORY)
    end subroutine new_gencorrs_fft

    subroutine new_order_corr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call order_corr%new(&
        &'order_corr',&                        ! name
        &'order_corr ',&                       ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call order_corr%add_input(UI_IO, )
        ! parameter input/output
        !call order_corr%add_input(UI_IMG, )
        ! <no additional inputs>
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
        call add_ui_program('order_corr', order_corr, tsttab, UI_CATEGORY)
    end subroutine new_order_corr

    subroutine new_phasecorr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call phasecorr%new(&
        &'phasecorr',&                         ! name
        &'phasecorr ',&                        ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call phasecorr%add_input(UI_IO, )
        ! parameter input/output
        !call phasecorr%add_input(UI_IMG, )
        ! <no additional inputs>
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
        call add_ui_program('phasecorr', phasecorr, tsttab, UI_CATEGORY)
    end subroutine new_phasecorr

    subroutine new_rank_weights( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call rank_weights%new(&
        &'rank_weights',&                      ! name
        &'rank_weights ',&                     ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call rank_weights%add_input(UI_IO, )
        ! parameter input/output
        !call rank_weights%add_input(UI_IMG, )
        ! <no additional inputs>
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
        call add_ui_program('rank_weights', rank_weights, tsttab, UI_CATEGORY)
    end subroutine new_rank_weights

    subroutine new_rotate_ref( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call rotate_ref%new(&
        &'rotate_ref',&                        ! name
        &'rotate_ref ',&                       ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call rotate_ref%add_input(UI_IO, )
        ! parameter input/output
        !call rotate_ref%add_input(UI_IMG, )
        ! <no additional inputs>
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
        call add_ui_program('rotate_ref', rotate_ref, tsttab, UI_CATEGORY)
    end subroutine new_rotate_ref

end module simple_test_ui_fft
