program simple_test_ui_navigation
use simple_string, only: string
use simple_test_utils, only: assert_char, assert_int, assert_true, report_summary, tests_failed
use simple_ui_navigation, only: UI_LAYOUT_FLAT, UI_LAYOUT_GROUPED, ui_layout_is_valid, ui_layout_name, &
                               &ui_program_group, ui_suite
implicit none

type(string)            :: program_names(2)
type(ui_program_group)  :: groups(2), empty_group
type(ui_suite)          :: suite

program_names(1) = 'refine3D'
program_names(2) = 'postprocess'
call groups(1)%new('refinement3d', '3D refinement', 60, program_names)
call groups(2)%new('validation', 'Validation', 70)
call suite%new('simple_exec', 'simple_exec', 'SIMPLE', UI_LAYOUT_GROUPED, groups)
call empty_group%new('workflows', 'Stream workflows', 10)

call assert_true(ui_layout_is_valid(UI_LAYOUT_GROUPED), 'grouped layout is valid')
call assert_true(ui_layout_is_valid(UI_LAYOUT_FLAT), 'flat layout is valid')
call assert_true(.not. ui_layout_is_valid(0), 'zero layout is invalid')
call assert_char('grouped', ui_layout_name(UI_LAYOUT_GROUPED), 'grouped layout name')
call assert_char('flat', ui_layout_name(UI_LAYOUT_FLAT), 'flat layout name')
call assert_char('invalid', ui_layout_name(0), 'invalid layout name')

call assert_char('simple_exec', suite%id%to_char(), 'suite identifier')
call assert_char('simple_exec', suite%executable%to_char(), 'suite executable')
call assert_char('SIMPLE', suite%display_name%to_char(), 'suite display name')
call assert_int(UI_LAYOUT_GROUPED, suite%layout, 'suite layout')
call assert_int(2, suite%ngroups(), 'suite group count')
call assert_char('refinement3d', suite%groups(1)%id%to_char(), 'group identifier')
call assert_char('3D refinement', suite%groups(1)%title%to_char(), 'group title')
call assert_int(60, suite%groups(1)%display_order, 'group display order')
call assert_int(2, suite%groups(1)%nprograms(), 'group program count')
call assert_char('postprocess', suite%groups(1)%program_names(2)%to_char(), 'group program name')
call assert_int(0, empty_group%nprograms(), 'empty group has no programs')

call report_summary()
if( tests_failed /= 0 ) error stop 1

end program simple_test_ui_navigation
