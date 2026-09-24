!@descr: exact-value unit tests for ANSI control-string formatting
module simple_ansi_ctrls_tester
use simple_ansi_ctrls, only: C_BLACK, C_RED, C_GREEN, C_YELLOW, C_BLUE, C_MAGENTA, C_CYAN, C_WHITE, &
    &C_MARKED_BLACK, C_MARKED_RED, C_MARKED_GREEN, C_MARKED_YELLOW, C_MARKED_BLUE, C_MARKED_MAGENTA, &
    &C_MARKED_CYAN, C_MARKED_WHITE, C_BOLD, format_str
use simple_test_utils, only: assert_char, assert_int
implicit none
private
public :: run_all_ansi_ctrls_tests

contains

    subroutine run_all_ansi_ctrls_tests
        character(len=1), parameter :: ESC = achar(27)
        character(len=:), allocatable :: formatted

        call assert_char('30', C_BLACK,   'ANSI black foreground code')
        call assert_char('31', C_RED,     'ANSI red foreground code')
        call assert_char('32', C_GREEN,   'ANSI green foreground code')
        call assert_char('33', C_YELLOW,  'ANSI yellow foreground code')
        call assert_char('34', C_BLUE,    'ANSI blue foreground code')
        call assert_char('35', C_MAGENTA, 'ANSI magenta foreground code')
        call assert_char('36', C_CYAN,    'ANSI cyan foreground code')
        call assert_char('37', C_WHITE,   'ANSI white foreground code')

        call assert_char('40', C_MARKED_BLACK,   'ANSI black background code')
        call assert_char('41', C_MARKED_RED,     'ANSI red background code')
        call assert_char('42', C_MARKED_GREEN,   'ANSI green background code')
        call assert_char('43', C_MARKED_YELLOW,  'ANSI yellow background code')
        call assert_char('44', C_MARKED_BLUE,    'ANSI blue background code')
        call assert_char('45', C_MARKED_MAGENTA, 'ANSI magenta background code')
        call assert_char('46', C_MARKED_CYAN,    'ANSI cyan background code')
        call assert_char('47', C_MARKED_WHITE,   'ANSI white background code')

        formatted = format_str('Red', C_RED)
        call assert_char(ESC//'[31mRed'//ESC//'[0m', formatted, &
            &'format_str emits red start code, unchanged text, and reset code')
        call assert_int(12, len(formatted), 'formatted red text has the exact ANSI sequence length')

        call assert_char(ESC//'[30mBlack'//ESC//'[0m', format_str('Black', C_BLACK), &
            &'format_str emits the exact black foreground sequence')
        call assert_char(ESC//'[32mGreen'//ESC//'[0m', format_str('Green', C_GREEN), &
            &'format_str emits the exact green foreground sequence')
        call assert_char(ESC//'[33mYellow'//ESC//'[0m', format_str('Yellow', C_YELLOW), &
            &'format_str emits the exact yellow foreground sequence')
        call assert_char(ESC//'[34mBlue'//ESC//'[0m', format_str('Blue', C_BLUE), &
            &'format_str emits the exact blue foreground sequence')
        call assert_char(ESC//'[35mMagenta'//ESC//'[0m', format_str('Magenta', C_MAGENTA), &
            &'format_str emits the exact magenta foreground sequence')
        call assert_char(ESC//'[36mCyan'//ESC//'[0m', format_str('Cyan', C_CYAN), &
            &'format_str emits the exact cyan foreground sequence')
        call assert_char(ESC//'[37mWhite'//ESC//'[0m', format_str('White', C_WHITE), &
            &'format_str emits the exact white foreground sequence')

        call assert_char(ESC//'[40mBlack'//ESC//'[0m', format_str('Black', C_MARKED_BLACK), &
            &'format_str emits the exact black background sequence')
        call assert_char(ESC//'[41mRed'//ESC//'[0m', format_str('Red', C_MARKED_RED), &
            &'format_str emits the exact red background sequence')
        call assert_char(ESC//'[42mGreen'//ESC//'[0m', format_str('Green', C_MARKED_GREEN), &
            &'format_str emits the exact green background sequence')
        call assert_char(ESC//'[43mYellow'//ESC//'[0m', format_str('Yellow', C_MARKED_YELLOW), &
            &'format_str emits the exact yellow background sequence')
        call assert_char(ESC//'[44mBlue'//ESC//'[0m', format_str('Blue', C_MARKED_BLUE), &
            &'format_str emits the exact blue background sequence')
        call assert_char(ESC//'[45mMagenta'//ESC//'[0m', format_str('Magenta', C_MARKED_MAGENTA), &
            &'format_str emits the exact magenta background sequence')
        call assert_char(ESC//'[46mCyan'//ESC//'[0m', format_str('Cyan', C_MARKED_CYAN), &
            &'format_str emits the exact cyan background sequence')
        call assert_char(ESC//'[47mWhite'//ESC//'[0m', format_str('White', C_MARKED_WHITE), &
            &'format_str emits the exact white background sequence')

        formatted = format_str(' Alert ', C_BOLD//';'//C_RED)
        call assert_char(ESC//'[1;31m Alert '//ESC//'[0m', formatted, &
            &'format_str accepts a combined bold-red code and preserves spaces')
        call assert_int(18, len(formatted), 'combined formatting has the exact ANSI sequence length')
    end subroutine run_all_ansi_ctrls_tests

end module simple_ansi_ctrls_tester
