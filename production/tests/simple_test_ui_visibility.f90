program simple_test_ui_visibility
use simple_test_utils,    only: assert_char, assert_int, assert_true, report_summary, tests_failed
use simple_linked_list,   only: linked_list, list_iterator
use simple_ui_param,      only: ui_param
use simple_ui_program,    only: UI_PARM, ui_program
use simple_ui_visibility, only: UI_VIS_STANDARD, UI_VIS_ADVANCED, UI_VIS_DEVELOPER, &
                               &ui_visibility_from_advanced, ui_visibility_is_valid, ui_visibility_name
implicit none

type(ui_param)   :: param
type(ui_program) :: ui_prg

call assert_int(UI_VIS_STANDARD, ui_visibility_from_advanced(.false.), 'legacy standard mapping')
call assert_int(UI_VIS_ADVANCED, ui_visibility_from_advanced(.true.),  'legacy advanced mapping')
call assert_true_all_valid
call assert_char('standard',  ui_visibility_name(UI_VIS_STANDARD),  'standard visibility name')
call assert_char('advanced',  ui_visibility_name(UI_VIS_ADVANCED),  'advanced visibility name')
call assert_char('developer', ui_visibility_name(UI_VIS_DEVELOPER), 'developer visibility name')

call param%set_param('test_param', 'num', 'Test parameter', 'Test parameter help', 'e.g. 1', .false., 1.)
call assert_int(UI_VIS_DEVELOPER, param%visibility, 'unspecified parameter visibility is developer')
call param%apply_gui_overrides(gui_advanced=.false.)
call assert_int(UI_VIS_STANDARD, param%visibility, 'legacy parameter standard visibility')
call param%apply_gui_overrides(gui_visibility=UI_VIS_DEVELOPER)
call assert_int(UI_VIS_DEVELOPER, param%visibility, 'explicit parameter developer visibility')
call assert_true(param%has_default, 'optional numeric parameter has a default')
call assert_int(0, size(param%choices), 'numeric parameter has no choices')

call param%set_param('choice_param', 'multi', 'Choice parameter', 'Choice parameter help', &
    &'(first|second|third){second}', .false., 'second')
call assert_true(param%has_default, 'optional multi parameter has a default')
call assert_int(3, size(param%choices), 'multi parameter choices are structured')
call assert_char('first',  param%choices(1)%value%to_char(), 'first choice value')
call assert_char('second', param%choices(2)%label%to_char(), 'choice label defaults to CLI value')

call param%set_param('binary_param', 'binary', 'Binary parameter', 'Binary parameter help', &
    &'(yes|no){no}', .true., '')
call assert_true(.not. param%has_default, 'required binary parameter has no default')
call assert_int(2, size(param%choices), 'binary parameter choices are structured')

call ui_prg%new('test_program', 'Test program', 'Test program help', 'simple_exec', .false., &
    &gui_advanced=.false.)
call assert_int(UI_VIS_STANDARD, ui_prg%visibility, 'legacy program standard visibility')
call ui_prg%new('test_program', 'Test program', 'Test program help', 'simple_exec', .false., &
    &gui_visibility=UI_VIS_DEVELOPER)
call assert_int(UI_VIS_DEVELOPER, ui_prg%visibility, 'explicit program developer visibility')
call ui_prg%add_input(UI_PARM, 'visibility_param', 'num', 'Visibility parameter', &
    &'Visibility parameter help', 'e.g. 1', .false., 1., gui_visibility=UI_VIS_STANDARD)
call assert_input_visibility(ui_prg%parm_ios, 'visibility_param', UI_VIS_STANDARD)
call ui_prg%add_input(UI_PARM, param, descr_placeholder_override='(on|off){on}')
call assert_input_choice(ui_prg%parm_ios, 'binary_param', 'off')

call report_summary()
if( tests_failed /= 0 ) error stop 1

contains

    subroutine assert_true_all_valid
        use simple_test_utils, only: assert_true
        call assert_true(ui_visibility_is_valid(UI_VIS_STANDARD),  'standard visibility is valid')
        call assert_true(ui_visibility_is_valid(UI_VIS_ADVANCED),  'advanced visibility is valid')
        call assert_true(ui_visibility_is_valid(UI_VIS_DEVELOPER), 'developer visibility is valid')
        call assert_true(.not. ui_visibility_is_valid(0),          'zero visibility is invalid')
    end subroutine assert_true_all_valid

    subroutine assert_input_visibility( inputs, key, expected_visibility )
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key
        integer,           intent(in) :: expected_visibility
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_param )
                    if( input%key%to_char() == key )then
                        call assert_int(expected_visibility, input%visibility, 'explicit input visibility')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
                class default
                    call assert_int(1, 0, 'UI input is a ui_param')
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input exists')
    end subroutine assert_input_visibility

    subroutine assert_input_choice( inputs, key, expected_choice )
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key, expected_choice
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_param )
                    if( input%key%to_char() == key )then
                        call assert_char(expected_choice, input%choices(2)%value%to_char(), &
                            &'placeholder override refreshes structured choices')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
                class default
                    call assert_int(1, 0, 'UI input is a ui_param')
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input exists')
    end subroutine assert_input_choice

end program simple_test_ui_visibility
