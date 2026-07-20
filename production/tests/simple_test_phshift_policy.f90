program simple_test_phshift_policy
use simple_string,      only: string
use simple_linked_list, only: linked_list, list_iterator
use simple_ui,          only: make_ui, get_prg_ptr
use simple_ui_param,    only: ui_param
use simple_ui_program,  only: ui_program
implicit none

character(len=16), parameter :: FITTING_PROGRAMS(5) = [character(len=16) :: &
    &'ctf_estimate', 'preprocess', 'preproc', 'mini_stream', 'check_refpick']
type(ui_program), pointer :: program_ptr
integer :: i

call make_ui
do i = 1,size(FITTING_PROGRAMS)
    call get_prg_ptr(string(FITTING_PROGRAMS(i)), program_ptr)
    if( .not. associated(program_ptr) )then
        write(*,'(A)') 'Missing UI program: '//trim(FITTING_PROGRAMS(i))
        error stop 1
    endif
    call assert_ui_param(program_ptr%parm_ios, 'fit_phshift', FITTING_PROGRAMS(i), &
        &expected_type='binary', expected_default='no')
    call assert_ui_param(program_ptr%srch_ctrls, 'phshift_min',  FITTING_PROGRAMS(i))
    call assert_ui_param(program_ptr%srch_ctrls, 'phshift_max',  FITTING_PROGRAMS(i))
    call assert_ui_param(program_ptr%srch_ctrls, 'phshift_step', FITTING_PROGRAMS(i))
enddo

write(*,'(A,I0,A)') ' ALL TESTS PASSED (', size(FITTING_PROGRAMS), ' fitting-program UI contracts).'

contains

    subroutine assert_ui_param(params, key, program_name, expected_type, expected_default)
        type(linked_list), intent(in) :: params
        character(len=*),  intent(in) :: key, program_name
        character(len=*),  intent(in), optional :: expected_type, expected_default
        type(list_iterator)   :: iterator
        class(*), allocatable :: value
        logical :: found
        found    = .false.
        iterator = params%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type(param => value)
                type is(ui_param)
                    if( param%key%to_char() == key )then
                        found = .true.
                        if( present(expected_type) )then
                            if( param%keytype%to_char() /= expected_type )then
                                write(*,'(A)') trim(program_name)//': '//key//' has wrong UI type'
                                error stop 1
                            endif
                        endif
                        if( present(expected_default) )then
                            if( param%cval_default%to_char() /= expected_default )then
                                write(*,'(A)') trim(program_name)//': '//key//' has wrong default'
                                error stop 1
                            endif
                        endif
                    endif
                class default
                    write(*,'(A)') trim(program_name)//': invalid UI parameter-list entry'
                    error stop 1
            end select
            if( allocated(value) ) deallocate(value)
            if( found ) exit
            call iterator%next()
        enddo
        if( .not. found )then
            write(*,'(A)') trim(program_name)//': missing UI parameter '//key
            error stop 1
        endif
    end subroutine assert_ui_param

end program simple_test_phshift_policy
