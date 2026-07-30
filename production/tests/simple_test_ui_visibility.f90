program simple_test_ui_visibility
use simple_test_utils,    only: assert_char, assert_int, assert_true, report_summary, tests_failed
use simple_linked_list,   only: linked_list, list_iterator
use simple_string,        only: string
use simple_ui,            only: make_ui, make_test_ui, get_prg_ptr, get_test_prg_ptr
use simple_ui_param,      only: UI_PLACEHOLDER_MAX_LEN, ui_param
use simple_ui_program,    only: UI_DISPLAY_NAME_MAX_LEN, UI_FILE, UI_PARM, UI_SUMMARY_MAX_LEN, &
    &category_descriptor, ui_program, ui_program_input, ui_activation_equals_any
use simple_ui_descriptor_types, only: ui_choices
use simple_ui_visibility, only: UI_VIS_STANDARD, UI_VIS_ADVANCED, UI_VIS_DEVELOPER, &
                               &ui_visibility_is_valid, ui_visibility_name
implicit none

type(ui_param)   :: param
type(ui_program) :: ui_prg
type(ui_program), pointer :: registered_prg => null()
type(string) :: program_name
type(string), allocatable :: supplied_keys(:)

call assert_true_all_valid
call assert_char('standard',  ui_visibility_name(UI_VIS_STANDARD),  'standard visibility name')
call assert_char('advanced',  ui_visibility_name(UI_VIS_ADVANCED),  'advanced visibility name')
call assert_char('developer', ui_visibility_name(UI_VIS_DEVELOPER), 'developer visibility name')

call param%set_param('test_param', 'num', 'Test parameter', 'Test parameter help', 'e.g. 1', .false., 1.)
call assert_char('Test parameter', param%label%to_char(), 'input short text is the label')
call assert_char('Test parameter help', param%help%to_char(), 'input long text is help')
call assert_char('e.g. 10', param%placeholder%to_char(), 'numeric placeholder is a standard example')
call assert_true(len_trim(param%placeholder%to_char()) <= UI_PLACEHOLDER_MAX_LEN, &
    &'placeholder respects the display limit')
call assert_true(param%has_default, 'optional numeric parameter has a default')
call assert_int(0, size(param%choices), 'numeric parameter has no choices')

call param%set_param('choice_param', 'multi', 'Choice parameter', 'Choice parameter help', &
    &'', .false., 'second', choices=ui_choices([character(len=6) :: 'first', 'second', 'third']))
call assert_true(param%has_default, 'optional multi parameter has a default')
call assert_int(3, size(param%choices), 'multi parameter choices are structured')
call assert_char('first',  param%choices(1)%value%to_char(), 'first choice value')
call assert_char('second', param%choices(2)%label%to_char(), 'choice label defaults to CLI value')
call assert_char('', param%placeholder%to_char(), 'choice fields have no placeholder')

call param%set_param('binary_param', 'binary', 'Binary parameter', 'Binary parameter help', &
    &'', .true., '', choices=ui_choices([character(len=3) :: 'yes', 'no']))
call assert_true(.not. param%has_default, 'required binary parameter has no default')
call assert_int(2, size(param%choices), 'binary parameter choices are structured')

call ui_prg%new('test_program', 'Test program', 'Test program help', 'simple_exec', .false., &
    &visibility=UI_VIS_STANDARD)
call assert_char('Test program', ui_prg%summary%to_char(), 'program short text is the summary')
call assert_true(len_trim(ui_prg%summary%to_char()) <= UI_SUMMARY_MAX_LEN, &
    &'program summary respects the display limit')
call assert_char('Test program', ui_prg%display_name%to_char(), 'program display name falls back to the summary')
call assert_true(len_trim(ui_prg%display_name%to_char()) <= UI_DISPLAY_NAME_MAX_LEN, &
    &'program display name respects the display limit')
call assert_char('Test program help', ui_prg%help%to_char(), 'program long text is help')
call assert_int(UI_VIS_STANDARD, ui_prg%visibility, 'explicit program standard visibility')
call ui_prg%set_category(category_descriptor('test', 'Test Programs', 1))
call assert_char('test', ui_prg%category%to_char(), 'program category')
call assert_char('Test Programs', ui_prg%category_display_name%to_char(), 'program category display name')
call assert_int(1, ui_prg%category_order, 'program category order')
call ui_prg%new('test_program', 'Test program', 'Test program help', 'simple_exec', .false., &
    &visibility=UI_VIS_DEVELOPER, display_name='Visible test program')
call assert_char('Visible test program', ui_prg%display_name%to_char(), 'explicit program display name')
call assert_int(UI_VIS_DEVELOPER, ui_prg%visibility, 'explicit program developer visibility')
call ui_prg%add_input(UI_PARM, 'visibility_param', 'num', 'Visibility parameter', &
    &'Visibility parameter help', 'e.g. 1', .false., 1., group='test controls', &
    &visibility=UI_VIS_STANDARD, activation=ui_activation_equals_any('mode', &
    &[character(len=8) :: 'standard', 'expert']))
call assert_input_visibility(ui_prg%parm_ios, 'visibility_param', UI_VIS_STANDARD)
call assert_input_binding(ui_prg%parm_ios, 'visibility_param')
call ui_prg%add_input(UI_PARM, param, choices_override=ui_choices([character(len=3) :: 'on', 'off']))
call assert_input_choice(ui_prg%parm_ios, 'binary_param', 'off')
call ui_prg%add_input(UI_FILE, 'input_a', 'file', 'Input A', 'First input source', 'e.g. input-a.mrc', .false., '')
call ui_prg%add_input(UI_FILE, 'input_b', 'file', 'Input B', 'Second input source', 'e.g. input-b.mrc', .false., '')
call assert_input_visibility(ui_prg%file_ios, 'input_a', UI_VIS_DEVELOPER)
call ui_prg%add_requirement('input_source', 'Input source', 'Supply exactly one input source.', &
    &[character(len=7) :: 'input_a', 'input_b'], max_selected=1)
call assert_int(1, size(ui_prg%requirements), 'one requirement group is registered')
call assert_int(1, ui_prg%requirements(1)%min_selected, 'requirement minimum selection count')
call assert_int(1, ui_prg%requirements(1)%max_selected, 'requirement maximum selection count')
allocate(supplied_keys(1))
supplied_keys(1) = 'input_a'
call assert_true(ui_prg%requirements_satisfied(supplied_keys), 'one accepted input satisfies requirement')
supplied_keys(1) = 'other'
call assert_true(.not. ui_prg%requirements_satisfied(supplied_keys), 'unrelated input does not satisfy requirement')
deallocate(supplied_keys)

call make_ui
call assert_registered_requirement('binarize', 'input', 1, 1)
call assert_registered_category('icm2D', 'denoise', 'Denoising', 70)
call assert_registered_category('atoms_stats', 'atom', 'Atom Analysis', 50)
call assert_registered_category('abinitio2D_stream', 'stream', 'Stream Workflows', 10)
call make_test_ui
program_name = 'strategy2D'
call get_test_prg_ptr(program_name, registered_prg)
call assert_true(associated(registered_prg), 'test program is registered')
if( associated(registered_prg) )then
    call assert_char('class', registered_prg%category%to_char(), 'test program category')
endif

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

    subroutine assert_registered_category( name, expected_category, expected_display_name, expected_order )
        character(len=*), intent(in) :: name, expected_category, expected_display_name
        integer,          intent(in) :: expected_order
        program_name = name
        call get_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), trim(name)//' is registered')
        if( associated(registered_prg) )then
            call assert_char(expected_category, registered_prg%category%to_char(), trim(name)//' category')
            call assert_char(expected_display_name, registered_prg%category_display_name%to_char(), &
                &trim(name)//' category display name')
            call assert_int(expected_order, registered_prg%category_order, trim(name)//' category order')
        endif
    end subroutine assert_registered_category

    subroutine assert_registered_requirement(name, expected_id, expected_minimum, expected_maximum)
        character(len=*), intent(in) :: name, expected_id
        integer, intent(in) :: expected_minimum, expected_maximum
        program_name = name
        call get_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), trim(name)//' is registered')
        if (associated(registered_prg)) then
            call assert_int(1, size(registered_prg%requirements), trim(name)//' requirement count')
            if (size(registered_prg%requirements) == 1) then
                call assert_char(expected_id, registered_prg%requirements(1)%id%to_char(), &
                    &trim(name)//' requirement id')
                call assert_int(expected_minimum, registered_prg%requirements(1)%min_selected, &
                    &trim(name)//' requirement minimum')
                call assert_int(expected_maximum, registered_prg%requirements(1)%max_selected, &
                    &trim(name)//' requirement maximum')
            endif
        endif
    end subroutine assert_registered_requirement

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
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        call assert_int(expected_visibility, input%visibility, 'explicit input visibility')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
                class default
                    call assert_int(1, 0, 'UI input is a ui_program_input')
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
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        call assert_char(expected_choice, input%param%choices(2)%value%to_char(), &
                            &'explicit choice override replaces structured choices')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
                class default
                    call assert_int(1, 0, 'UI input is a ui_program_input')
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input exists')
    end subroutine assert_input_choice

    subroutine assert_input_binding( inputs, key )
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        call assert_char('test controls', input%group%id%to_char(), 'input group id')
                        call assert_int(1, input%group%order, 'first input group order')
                        call assert_char('mode', input%activation%key%to_char(), 'activation key')
                        call assert_int(2, size(input%activation%equals_any), 'activation value count')
                        call assert_char('expert', input%activation%equals_any(2)%to_char(), 'activation value')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input binding exists')
    end subroutine assert_input_binding

end program simple_test_ui_visibility
