!@descr: module defining the user interfaces for fft testprograms in the simple_test_exec suite
module simple_test_ui_fft
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('fft', 'FFT', 20)
type(ui_program), target :: gencorrs_fft

contains

    subroutine construct_test_fft_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_gencorrs_fft(tsttab)
    end subroutine construct_test_fft_programs

    subroutine new_gencorrs_fft( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call gencorrs_fft%new(&
        &'gencorrs_fft',&                      ! name
        &'polar correlation of generated images',& ! summary
        &'polarises generated images into a polarft_calc and checks gen_objfun_vals peaks (self, rotated copy, unrelated)',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! <no inputs: the test generates its own images>
        ! add to ui_hash
        call add_ui_program('gencorrs_fft', gencorrs_fft, tsttab, UI_CATEGORY)
    end subroutine new_gencorrs_fft

end module simple_test_ui_fft
