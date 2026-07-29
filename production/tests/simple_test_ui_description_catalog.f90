program simple_test_ui_description_catalog
use simple_test_utils, only: assert_char, assert_int, assert_true, report_summary, tests_failed
use simple_linked_list, only: linked_list, list_iterator
use simple_ui_descriptor_types, only: ui_choice
use simple_ui_description_catalog, only: get_ui_program_presentation, get_ui_input_presentation, &
                                         &apply_ui_input_choice_presentation
use simple_ui_param, only: ui_param
use simple_ui_program, only: UI_PARM, ui_program
use simple_ui_visibility, only: UI_VIS_STANDARD, UI_VIS_DEVELOPER
implicit none

character(len=:), allocatable :: display_name, summary, help, label, placeholder, units
integer :: visibility
logical :: found
type(ui_choice) :: choices(2)
type(ui_program) :: program

call get_ui_program_presentation('simple_exec', 'refine3D', display_name, summary, help, visibility, found)
call assert_true(found, 'refine3D program record exists in catalog')
call assert_char('3D refinement', display_name, 'catalog program display name')
call assert_char('3D refinement', summary, 'catalog program summary')
call assert_int(UI_VIS_STANDARD, visibility, 'catalog program visibility')

call get_ui_input_presentation('simple_exec', 'refine3D', 'pgrp', label, help, placeholder, units, visibility, found)
call assert_true(found, 'refine3D pgrp record exists in catalog')
call assert_char('Point-group symmetry', label, 'catalog input label')
call assert_int(UI_VIS_STANDARD, visibility, 'catalog input visibility')

call program%new('refine3D', 'legacy summary', 'legacy help', 'simple_exec', .true.)
call assert_char('3D refinement', program%display_name%to_char(), 'catalog display name is applied during program construction')
call assert_char('3D refinement', program%descr_short%to_char(), 'catalog is applied during program construction')
call assert_int(UI_VIS_STANDARD, program%visibility, 'catalog program visibility is applied')
call program%add_input(UI_PARM, 'pgrp', 'multi', 'legacy label', 'legacy help', '(c1|d1){c1}', .false., 'c1')
call assert_catalog_input(program%parm_ios, 'pgrp')

choices(1)%value = 'yes'
choices(2)%value = 'no'
call apply_ui_input_choice_presentation('simple_exec', 'refine3D', 'center', choices, found)
call assert_true(found, 'center choices exist in catalog')
call assert_char('yes', choices(1)%label%to_char(), 'catalog choice label')

call get_ui_program_presentation('simple_exec', 'not_a_program', display_name, summary, help, visibility, found)
call assert_true(.not. found, 'missing program has no catalog record')
call assert_int(UI_VIS_DEVELOPER, visibility, 'missing program is developer by default')

call report_summary()
if( tests_failed /= 0 ) error stop 1

contains

    subroutine assert_catalog_input(inputs, key)
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type(input => value)
                type is(ui_param)
                    if(input%key%to_char() == key)then
                        call assert_char('Point-group symmetry', input%descr_short%to_char(), &
                            &'catalog is applied during input construction')
                        call assert_int(UI_VIS_STANDARD, input%visibility, 'catalog input visibility is applied')
                        if(allocated(value)) deallocate(value)
                        return
                    endif
            end select
            if(allocated(value)) deallocate(value)
            call iterator%next()
        enddo
        call assert_true(.false., 'catalog input exists')
    end subroutine assert_catalog_input

end program simple_test_ui_description_catalog
