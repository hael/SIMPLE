module continuous_inplane_refine3D_policy_test
implicit none
private
public :: run_policy_gate

contains

subroutine run_policy_gate()
use simple_cmdline, only: cmdline
use simple_parameters, only: parameters
use simple_refine3D_strategy, only: strip_refine3D_search_only_args
use simple_ui, only: make_ui
implicit none

type(cmdline)    :: cline
type(cmdline)    :: child_cline
type(parameters) :: defaults

if( trim(defaults%inpl_cont) /= 'yes' ) &
    &error stop 'refine3D continuous in-plane default is not yes'
call make_ui
if( .not. has_search_input('refine3D', 'inpl_cont', expected_default='yes') ) &
    &error stop 'refine3D UI does not expose inpl_cont=yes'
if( .not. has_search_input('refine3D_auto', 'inpl_cont', expected_default='yes') ) &
    &error stop 'refine3D_auto UI does not expose inpl_cont=yes'
if( .not. has_search_input('refine3D_states', 'inpl_cont', expected_default='yes') ) &
    &error stop 'refine3D_states UI does not expose inpl_cont=yes'

call cline%set('prg', 'refine3D')
call cline%set('inpl_cont', 'yes')
child_cline = cline
call strip_refine3D_search_only_args(child_cline)
if( child_cline%defined('inpl_cont') ) &
    &error stop 'refine3D matcher-only option leaked into a child workflow'
call child_cline%kill
call cline%kill

write(*,'(a)') 'REFINE3D_INPL_CONT_POLICY DEFAULT/UI/CHILD_STRIPPING: PASS'
end subroutine run_policy_gate

logical function has_search_input(program_name, key, expected_default) result(found)
use simple_linked_list, only: list_iterator
use simple_string, only: string
use simple_ui, only: get_prg_ptr
use simple_ui_program, only: ui_program, ui_program_input
implicit none

character(len=*), intent(in) :: program_name, key
character(len=*), optional, intent(in) :: expected_default
type(ui_program), pointer :: program => null()
type(list_iterator) :: iterator
type(string) :: name
class(*), allocatable :: value

found = .false.
name = program_name
call get_prg_ptr(name, program)
call name%kill
if( .not. associated(program) ) return
iterator = program%srch_ctrls%begin()
do while( iterator%has_value() )
    call iterator%getter(value)
    select type(input => value)
    type is(ui_program_input)
        if( input%param%key%to_char() == key )then
            found = .true.
            if( present(expected_default) )then
                found = input%param%has_default .and. &
                    &input%param%cval_default%to_char() == expected_default
            endif
            if( allocated(value) ) deallocate(value)
            return
        endif
    end select
    if( allocated(value) ) deallocate(value)
    call iterator%next()
enddo
end function has_search_input

end module continuous_inplane_refine3D_policy_test
