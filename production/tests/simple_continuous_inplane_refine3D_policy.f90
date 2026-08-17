module continuous_inplane_refine3D_policy_test
implicit none
private
public :: run_policy_gate, run_policy_rejection

contains

subroutine run_policy_gate()
use simple_cmdline, only: cmdline
use simple_commanders_refine3D, only: validate_refine3D_inpl_cont
use simple_defs, only: STDLEN
use simple_parameters, only: parameters, refine3D_inpl_cont_policy_error
use simple_refine3D_strategy, only: strip_refine3D_search_only_args
use simple_ui, only: make_ui
implicit none

type(cmdline)    :: cline
type(cmdline)    :: child_cline
type(parameters) :: defaults
character(len=STDLEN) :: policy_error

if( trim(defaults%inpl_cont) /= 'yes' ) &
    &error stop 'refine3D continuous in-plane default is not yes'
call make_ui
if( .not. has_search_input('refine3D', 'inpl_cont', expected_default='yes') ) &
    &error stop 'refine3D UI does not expose inpl_cont=yes'
if( .not. has_search_input('refine3D_auto', 'inpl_cont', expected_default='yes') ) &
    &error stop 'refine3D_auto UI does not expose inpl_cont=yes'
if( .not. has_search_input('refine3D_multi', 'inpl_cont', expected_default='yes') ) &
    &error stop 'refine3D_multi UI does not expose inpl_cont=yes'

policy_error = refine3D_inpl_cont_policy_error('no', 'cc', 'yes', &
    &'den', 'yes')
if( len_trim(policy_error) > 0 ) &
    &error stop 'shared refine3D policy altered an explicit default-off route'
policy_error = refine3D_inpl_cont_policy_error('yes', 'euclid', 'no', &
    &'raw', 'no')
if( len_trim(policy_error) > 0 ) &
    &error stop 'shared refine3D policy rejected the eligible opt-in route'
policy_error = refine3D_inpl_cont_policy_error('yes', 'cc', 'no', &
    &'raw', 'no')
if( len_trim(policy_error) > 0 ) &
    &error stop 'shared refine3D policy rejected the eligible cc opt-in route'
call cline%set('prg', 'refine3D')
call validate_refine3D_inpl_cont(cline)
call cline%set('inpl_cont', 'no')
call validate_refine3D_inpl_cont(cline)
call cline%set('inpl_cont', 'yes')
call cline%set('objfun', 'euclid')
call cline%set('objfun_den', 'no')
call cline%set('ptcl_src', 'raw')
call cline%set('projrec', 'no')
call cline%set('refine', 'shc')
call validate_refine3D_inpl_cont(cline)
call cline%set('refine', 'prob_neigh')
call validate_refine3D_inpl_cont(cline)
child_cline = cline
call strip_refine3D_search_only_args(child_cline)
if( child_cline%defined('inpl_cont') ) &
    &error stop 'refine3D matcher-only option leaked into a child workflow'
call child_cline%kill
call cline%kill

write(*,'(a)') 'REFINE3D_INPL_CONT_POLICY SHARED/DEFAULT/OFF/ELIGIBLE: PASS'
end subroutine run_policy_gate

subroutine run_policy_rejection(label)
use simple_cmdline, only: cmdline
use simple_commanders_refine3D, only: validate_refine3D_inpl_cont
implicit none

character(len=*), intent(in) :: label
type(cmdline) :: cline

call cline%set('prg', 'refine3D')
call cline%set('inpl_cont', 'yes')
call cline%set('objfun', 'euclid')
call cline%set('objfun_den', 'no')
call cline%set('ptcl_src', 'raw')
call cline%set('projrec', 'no')
call cline%set('refine', 'shc')
select case(trim(label))
case('policy_bad_value')
    call cline%set('inpl_cont', 'invalid')
case('policy_hybrid')
    call cline%set('objfun_den', 'yes')
case('policy_denoised')
    call cline%set('ptcl_src', 'den')
case('policy_projrec')
    call cline%set('projrec', 'yes')
case default
    error stop 'unknown refine3D policy rejection case'
end select
call validate_refine3D_inpl_cont(cline)
error stop 'unsupported refine3D continuous in-plane command was accepted'
end subroutine run_policy_rejection

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
