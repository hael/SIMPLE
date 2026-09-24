!@descr: module defining the user interface of the stream workflow test program (preproc) in the simple_test_exec suite
module simple_test_ui_stream
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('stream', 'Stream', 130)
type(ui_program), target :: preproc

contains

    subroutine construct_test_stream_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_preproc(tsttab)
    end subroutine construct_test_stream_programs

    subroutine new_preproc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call preproc%new('preproc', 'Streaming preprocessing', &
            &'generates five synthetic movies and validates streaming motion-correction and CTF outputs', &
            &'simple_test_exec', .false.)
        call add_ui_program('preproc', preproc, tsttab, UI_CATEGORY)
    end subroutine new_preproc

end module simple_test_ui_stream
