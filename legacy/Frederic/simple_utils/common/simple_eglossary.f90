! ============================================================================
! Name        : simple_glossary
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 27th of June 2016
! Description : Module to define a glossary of fucntionalities, thisd includes
!             : the simple_error_handling module, the simple_yaml module and 
!             : the like. This  module is a higher level on the
!             : simple_error_handling module which handles the error separately
! ============================================================================
!
module simple_eglossary
  use simple_err_defs
  use simple_error_handling
  use simple_eglossary_lowlev
  use simple_yaml_strings

  implicit none
!  character(len=max_field_length), parameter :: NOT_A_VALUE="__not_a_value__"

  private 

  type, public :: list_container
     character(len=max_field_length) :: val=' '
     type(eglossary),pointer :: egloss=>null() 
  end type list_container
  
  type, public :: eglossary_container
     character(len=max_field_length) :: key=' '
     character(len=max_field_length) :: value=' '
     type(eglossary), pointer :: child=>null()
  end type eglossary_container

  !Error codes
  integer, save, public :: EGLOSS_VALUE_ABSENT
  integer, save, public :: EGLOSS_INVALID_LIST
  integer, save, public :: EGLOSS_INVALID
  integer, save, public :: EGLOSS_CONVERSION_ERROR

  !control the error environment in simple_error_handling
  type(eglossary),pointer :: egloss_present_error=>null()

  interface assignment(=)
     !TODO: implement the interface with the methods and add as required
     module procedure get_value, get_integer, get_real, get_double
  end interface assignment(=)

  interface operator(.get.)
     module procedure list_container_if_key_exists
  end interface operator(.get.)

  interface operator(.in.)
     module procedure key_in_eglossary
  end interface operator(.in.)

  interface operator(.is.)
     !TODO: need to implement the is interface
     module procedure egloss_cont_new_with_value,egloss_cont_new_with_int
     module procedure egloss_cont_new_with_real,egloss_cont_new_with_double
     module procedure egloss_cont_new_with_egloss
     module procedure egloss_cont_new_with_double_v,egloss_cont_new_with_real_v
     module procedure egloss_cont_new_with_int_v, egloss_cont_new_with_value_v
  end interface operator(.is.)

  interface operator(.item.)
     module procedure item_char,item_double,item_real,item_integer,item_egloss
  end interface operator(.item.)

  interface list_new
     module procedure list_new, list_new_elems
  end interface list_new
  
  interface operator(.index.)
     module procedure find_index
  end interface operator(.index.)
  
  interface egloss_remove
     module procedure remove_egloss, remove_item
  end interface egloss_remove

  interface egloss_new
     module procedure egloss_new,egloss_new_elems
  end interface egloss_new
  
  interface egloss_iter
     module procedure egloss_iter,egloss_iter_lc
  end interface egloss_iter

  interface set
     module procedure put_value, put_child, put_list 
     module procedure put_integer,put_real, put_double, put_long
     module procedure put_logical
     module procedure put_listd,put_listi
  end interface set

  interface add
     !TODO: need to impolement the handlers for the interafce
     !module procedure add_real,charcater,integer, list , etc...
     module procedure add_char, add_egloss,add_integer,add_real,add_double
     module procedure add_long,add_log
  end interface add
     
  !overloaded operators
  public :: assignment(=)
  public :: operator(//),operator(.is.)
  public :: operator(.get.),operator(.in.),operator(.index.)
  !intertface for main handlers
  public :: set,add
  !Construcor
  public :: egloss_new
  !low level methods
  public :: egloss_init
  !handlers
  public :: has_key,glossary_errors,egloss_appender
  !public variables
  public :: eglossary
  !interfaced methods
  public :: file_get_last_error
  !iterators
  public :: egloss_iter,egloss_next
contains
  
  !glossary errors suites calls
  subroutine glossary_errors()
    implicit none
    !Sucess
    call file_err_define('SUCCESS', 'Operation has succeeded',&
         ERR_SUCCESS, err_action='No action')
    !Failure
    call file_err_define('FAIL', 'Operation has failed',&
         ERR_FAIL, err_action='stop')
    !Generic
    call file_err_define('GENERIC', 'Operation has not succeeded',&
         ERR_GENERIC, err_action='No action')
    !Undefined
    call file_err_define('ERR_NOT_DEFINED','The error id or name is invalid',&
         ERR_NOT_DEFINED, err_action='Check if the err id exists')

    return
  end subroutine glossary_errors
  !Pseudo constructor of the eglossary object
  function egloss_new(eglosses)
    type(eglossary_container), dimension(:), intent(in) :: eglosses
    type(eglossary), pointer :: egloss_new
    !local variables
    integer :: i_st,n_st
    type(eglossary), pointer :: egloss_tmp
    
    call egloss_init(egloss_tmp)
    n_st = size(eglosses)
    do i_st=1,n_st
       if (associated(eglosses(i_st)%child)) then
          !TODO: insert a setter for the eglosses and operator//
          !//eglosses(i_st)%key
          call set(egloss_tmp//eglosses(i_st)%key, eglosses(i_st)%child)
       else
          !TODO: insert a setter for the eglosses and operator//
          !//eglosses(i_st)%key
          call set(egloss_tmp//eglosses(i_st)%key, eglosses(i_st)%value)
       end if
          
    end do
    egloss_new => egloss_tmp
  end function egloss_new

  !Defines a eglossary from a array of storage data
  function egloss_new_elems(egloss0, egloss1, egloss2, egloss3, egloss4, &
                            egloss5, egloss6, egloss7, egloss8, egloss9, &
                            egloss10, egloss11, egloss12, egloss13, egloss14, &
                            egloss15, egloss16, egloss17, egloss18, egloss19)
    type(eglossary_container), intent(in), optional :: egloss0,egloss1,egloss2
    type(eglossary_container), intent(in), optional :: egloss3,egloss4,egloss5
    type(eglossary_container), intent(in), optional :: egloss6,egloss7,egloss8
    type(eglossary_container), intent(in), optional :: egloss9,egloss10,egloss11
    type(eglossary_container), intent(in), optional :: egloss12,egloss13
    type(eglossary_container), intent(in), optional :: egloss14,egloss15
    type(eglossary_container), intent(in), optional :: egloss16,egloss17
    type(eglossary_container), intent(in), optional :: egloss18,egloss19
    type(eglossary), pointer :: egloss_new_elems
    !local variables
    type(eglossary), pointer :: egloss_tmp

    call egloss_init(egloss_tmp)
    if (present(egloss0 )) call add_elem(egloss_tmp, egloss0 )
    if (present(egloss1 )) call add_elem(egloss_tmp, egloss1 )
    if (present(egloss2 )) call add_elem(egloss_tmp, egloss2 )
    if (present(egloss3 )) call add_elem(egloss_tmp, egloss3 )
    if (present(egloss4 )) call add_elem(egloss_tmp, egloss4 )
    if (present(egloss5 )) call add_elem(egloss_tmp, egloss5 )
    if (present(egloss6 )) call add_elem(egloss_tmp, egloss6 )
    if (present(egloss7 )) call add_elem(egloss_tmp, egloss7 )
    if (present(egloss8 )) call add_elem(egloss_tmp, egloss8 )
    if (present(egloss9 )) call add_elem(egloss_tmp, egloss9 )
    if (present(egloss10)) call add_elem(egloss_tmp, egloss10)
    if (present(egloss11)) call add_elem(egloss_tmp, egloss11)
    if (present(egloss12)) call add_elem(egloss_tmp, egloss12)
    if (present(egloss13)) call add_elem(egloss_tmp, egloss13)
    if (present(egloss14)) call add_elem(egloss_tmp, egloss14)
    if (present(egloss15)) call add_elem(egloss_tmp, egloss15)
    if (present(egloss16)) call add_elem(egloss_tmp, egloss16)
    if (present(egloss17)) call add_elem(egloss_tmp, egloss17)
    if (present(egloss18)) call add_elem(egloss_tmp, egloss18)
    if (present(egloss19)) call add_elem(egloss_tmp, egloss19)
    egloss_new_elems => egloss_tmp
  contains
    subroutine add_elem(egloss, elem)
      implicit none
      type(eglossary_container), intent(in) :: elem
      type(eglossary), pointer :: egloss
      
      if (associated(elem%child)) then
         call set(egloss//elem%key, elem%child)
      else
         call set(egloss//elem%key, elem%value)
      end if
    end subroutine add_elem
  end function egloss_new_elems

  !handlers
  !to append a eglossary to another
  recursive subroutine egloss_appender(egloss,brother)
    implicit none
    type(eglossary), pointer :: egloss
    type(eglossary), pointer :: brother
    !local variables
    type(eglossary), pointer :: egloss_tmp,iter

    if (.not. associated(brother)) return

    if (.not. associated(egloss)) then
       if (associated(brother%parent)) then
          call egloss_init(egloss)
          call set(egloss,brother)
       else
          egloss=>brother
       end if
    else if (.not. associated(egloss%parent)) then
       call egloss_appender(egloss%child,brother)
    else if (.not. associated(brother%parent)) then
       !increment the number of elements
       egloss%parent%data%nelems=&
            egloss%parent%data%nelems+brother%data%nelems
       call define_parent(egloss%parent,brother%child)
       call egloss_appender(egloss,brother%child)
       nullify(brother%child)
       call egloss_free(brother)
    else if (associated(egloss%previous)) then
       call egloss_appender(egloss%previous,brother)
    else
       egloss_tmp=>brother
       call egloss_appender(brother,egloss)
       egloss=>egloss_tmp
    end if
    
  end subroutine egloss_appender
  
  function has_key(egloss,key)
    !use simple_eglossary_lowlev
    implicit none
    type(eglossary), intent(in), pointer :: egloss
    character(len=*), intent(in) :: key
    logical :: has_key
    if (.not. associated(egloss) .or. trim(key)=="" )then
       has_key=.false.
       return
    end if
    has_key=associated(egloss%child)
    if (has_key) has_key=associated(get_egloss_from_the_key(egloss%child,key))
  end function has_key

  !puts value into eglossary
  subroutine put_value(egloss,val)
    use simple_eglossary_lowlev
    implicit none
    
    type(eglossary), pointer :: egloss
    character(len=*), intent(in) :: val
    !TODO: implement the insertion of the setter
    if (trim(val) == NOT_A_VALUE) then
       call file_err_throw('Invalid assignment for key "'//&
            trim(egloss%data%key)//'"',err_id=EGLOSS_VALUE_ABSENT)
       return
    end if
    write(*,*) "val: ",val
    write(*,*) "egloss%data%value: ",egloss%data%value
    
    if (associated(egloss%child) ) then
       call free_child(egloss)
    end if
    call file_strcpy(src=val,dest=egloss%data%value)
    return
  end subroutine put_value

  !true type insertors
  !integer
  subroutine put_integer(egloss,ival,fmt)
    use simple_yaml_strings
    implicit none
    type(eglossary),pointer :: egloss
    integer(kind=4),intent(in) :: ival
    character(len=*),optional, intent(in) :: fmt
    !TODO: need to have a put_value for real, int and double method
    write(*,*) "in put_integer"
    if (present(fmt)) then
       call put_value(egloss,trim(  adjustl( yaml_toa(ival,fmt=fmt) ) ) )
    else
       call put_value(egloss,adjustl( trim ( yaml_toa(ival) ) ) )
    end if
    return
  end subroutine put_integer
  !real
  subroutine put_real(egloss,rval,fmt)
    use simple_yaml_strings
    implicit none
    type(eglossary),pointer :: egloss
    real(kind=4),intent(in) :: rval
    character(len=*),optional, intent(in) :: fmt
    !TODO: need to have a put_value for real, int and double method
    write(*,*) "in put_real"
    if (present(fmt)) then
       call put_value(egloss,trim(  adjustl( yaml_toa(rval,fmt=fmt) ) ) )
    else
       call put_value(egloss,adjustl( trim ( yaml_toa(rval) ) ) )
    end if
    return
  end subroutine put_real
  
  !double
  subroutine put_double(egloss,dval,fmt)
    use simple_yaml_strings
    implicit none
    type(eglossary),pointer :: egloss
    real(kind=8),intent(in) :: dval
    character(len=*),optional, intent(in) :: fmt
    !TODO: need to have a put_value for real, int and double method
    write(*,*) "in put_double"
    if (present(fmt)) then
       call put_value(egloss,trim(  adjustl( yaml_toa(dval,fmt=fmt) ) ) )
    else
       call put_value(egloss,adjustl( trim ( yaml_toa(dval) ) ) )
    end if
    return
  end subroutine put_double
  !long
  subroutine put_long(egloss,lval,fmt)
    use simple_yaml_strings
    implicit none
    type(eglossary),pointer :: egloss
    integer(kind=8),intent(in) :: lval
    character(len=*),optional, intent(in) :: fmt
    !TODO: need to have a put_value for real, int and double method
    write(*,*) "in put_long"
    if (present(fmt)) then
       call put_value(egloss,trim(  adjustl( yaml_toa(lval,fmt=fmt) ) ) )
    else
       call put_value(egloss,adjustl( trim ( yaml_toa(lval) ) ) )
    end if
    return
  end subroutine put_long
  !TODO: need to implemewnt the logical version of this method writter
  !logical
  subroutine put_logical(egloss,tval,fmt)
    use simple_yaml_strings
    implicit none
    type(eglossary),pointer :: egloss
    logical,intent(in) :: tval
    character(len=*),optional, intent(in) :: fmt
    !TODO: need to have a put_value for real, int and double method
    write(*,*) "in put_logical"
    if (present(fmt)) then
       write(*,*) "TODO: need to iomplement the logical version of this"
       stop
!       call put_value(egloss,trim(  adjustl( yaml_toa(tval,fmt=fmt) ) ) )
    else
!       call put_value(egloss,adjustl( trim ( yaml_toa(tval) ) ) )
       write(*,*) "TODO: need to iomplement the logical version of this"
       stop
    end if
    return
  end subroutine put_logical
  
  !freeing the child node and setting the items en elems to 0
  subroutine free_child(egloss)
    implicit none
    type(eglossary), pointer :: egloss
    call egloss_free(egloss%child)
    !resetting the number of items
    egloss%data%nitems=0
    egloss%data%nelems=0
    return
  end subroutine free_child
  !puuting the 
  subroutine put_child(egloss,subd)
    use simple_eglossary_lowlev
    implicit none
    type(eglossary),pointer :: egloss
    type(eglossary),pointer :: subd
    !if the eglossary starts with a master tree, replace with child
    if (.not. associated(subd%parent) .and. associated(subd%child)) then
       call put_child(egloss,subd%child)
       nullify(subd%child)
       call egloss_free(subd)
       return
    end if

    call file_strcpy(src=' ',dest=egloss%data%value)
    if ( .not. associated(egloss%child,target=subd) .and. &
         associated(egloss%child)) then
       call free_child(egloss)
    end if
    egloss%child=>subd
    if (associated(subd%parent)) then
       !inherit the number of elements or items from subd's parent
       !which is guaranteed to be associated
       egloss%data%nelems=subd%parent%data%nelems
       egloss%data%nitems=subd%parent%data%nitems
    end if
    !TODO: need to fix the undefined access
    call define_parent(egloss,egloss%child)

    return
  end subroutine put_child

  subroutine put_list(egloss,list)
    implicit none
    character(len=*), dimension(:), intent(in) :: list
    type(eglossary), pointer :: egloss
    !local variables
    integer :: item,nitems,nitems_old

    nitems=size(list)
    nitems_old=egloss_len(egloss)
    do item=1,nitems
       !TODO: need to fix the setter using the overloaded operator
       !call set(egloss//(item-1),list(item))
    end do
    do item=nitems_old-1,nitems,-1
       call egloss_remove(egloss,item)
    end do

  end subroutine put_list

  subroutine put_listd(egloss,list)
    implicit none
    double precision, dimension(:), intent(in) :: list
    type(eglossary), pointer :: egloss
    !local variables
    integer :: item,nitems,nitems_old

    nitems=size(list)
    nitems_old=egloss_len(egloss)
    do item=1,nitems
       !call set(egloss//(item-1),list(item))
    end do
    do item=nitems_old-1,nitems,-1
       call egloss_remove(egloss,item)
    end do


  end subroutine put_listd

  subroutine put_listi(egloss,list)
    implicit none
    integer, dimension(:), intent(in) :: list
    type(eglossary), pointer :: egloss
    !local variables
    integer :: item,nitems,nitems_old

    nitems=size(list)
    nitems_old=egloss_len(egloss)
    do item=1,nitems
       !call set(egloss//(item-1),list(item))
    end do
    do item=nitems_old-1,nitems,-1
       !TODO: need to implement the interface for
       call egloss_remove(egloss,item)
    end do
  end subroutine put_listi

  
  subroutine remove_egloss(egloss,key,destroy)
    implicit none
    type(eglossary), pointer :: egloss 
    character(len=*), intent(in) :: key
    logical, intent(in), optional :: destroy
    !local variables

    !TODO: implement the method

  contains
    subroutine pop_egloss_(egloss,key,dst)
      implicit none
      type(eglossary), intent(inout), pointer :: egloss 
      character(len=*), intent(in) :: key
      logical, intent(in) :: dst
      !local variables
    end subroutine pop_egloss_

  end subroutine remove_egloss

  subroutine remove_item(egloss,item,destroy)
    implicit none
    type(eglossary), pointer :: egloss 
    integer, intent(in) :: item
    logical, intent(in), optional :: destroy
    !TODO: implement the method

  contains
    subroutine pop_item_(egloss,item,dst)
      implicit none
      type(eglossary), intent(inout), pointer :: egloss
      integer, intent(in) :: item
      logical, intent(in) :: dst
      !local variables
    end subroutine pop_item_

  end subroutine remove_item

  !implementation of the interface .get. operator
  function list_container_if_key_exists(egloss,key) result(list) 
    implicit none
    type(eglossary), pointer, intent(in) :: egloss
    character(len=*),intent(in) :: key
    type(list_container) :: list
    if(trim(key) .in. egloss ) list%egloss=>egloss//trim(key)
  end function list_container_if_key_exists

  !implementation of the .in. operator
  function key_in_eglossary(key,egloss)
    implicit none
    type(eglossary), pointer, intent(in) :: egloss
    character(len=*),intent(in) :: key
    logical :: key_in_eglossary
    if (egloss_len(egloss)>0)then
       key_in_eglossary = (egloss .index. key)>=0
    else
       key_in_eglossary = has_key(egloss,key)
    end if
  end function key_in_eglossary

  !implementation of the .index. operator
  function find_index(egloss,name)
    implicit none
    type(eglossary),pointer,intent(in) :: egloss
    character(len=*),intent(in) :: name
    integer :: find_index
    !local variables
    integer :: ind
    type(eglossary),pointer :: egloss_tmp
    find_index=-1
    ind=-1
    if (associated(egloss)) then
       egloss_tmp=>egloss_iter(egloss)
       loop_find: do while(associated(egloss_tmp))
          ind=ind+1
          if (name_is(egloss_tmp,name)) then
             find_index=ind
             exit loop_find
          end if
          egloss_tmp=>egloss_next(egloss_tmp)
       end do loop_find
    end if
  end function find_index

  !implementation of the .is. operator
  !character
  function egloss_cont_new_with_value(key,val) result(cont)
    implicit none
    character(len=*), intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont%key(1:max_field_length) = key
    cont%value(1:max_field_length) = yaml_toa(val)
  end function egloss_cont_new_with_value
  !integer
  function egloss_cont_new_with_int(key,val) result(cont)
    implicit none
    integer, intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont%key(1:max_field_length) = key
    cont%value(1:max_field_length) = yaml_toa(val)
  end function egloss_cont_new_with_int
  !real
  function egloss_cont_new_with_real(key,val) result(cont)
    implicit none
    real, intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont%key(1:max_field_length) = key
    cont%value(1:max_field_length) = yaml_toa(val)
  end function egloss_cont_new_with_real
  !double
  function egloss_cont_new_with_double(key,val) result(cont)
    implicit none
    double precision, intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont%key(1:max_field_length) = key
    cont%value(1:max_field_length) = yaml_toa(val)
  end function egloss_cont_new_with_double
  !eglossary
  function egloss_cont_new_with_egloss(key, val)
    implicit none
    character(len = *), intent(in) :: key
    type(eglossary), pointer, intent(in) :: val
    type(eglossary_container) :: egloss_cont_new_with_egloss
    egloss_cont_new_with_egloss%key(1:max_field_length) = key
    egloss_cont_new_with_egloss%child => val
  end function egloss_cont_new_with_egloss
  !vectors of characters, double, integers and real
  function egloss_cont_new_with_value_v(key, val) result(cont)
    implicit none
    character(len = *), dimension(:), intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont=egloss_cont_new_with_egloss(key,list_new(.item. val))
  end function egloss_cont_new_with_value_v

  function egloss_cont_new_with_int_v(key, val) result(cont)
    implicit none
    integer, dimension(:), intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont=egloss_cont_new_with_egloss(key,list_new(.item. val))
  end function egloss_cont_new_with_int_v

  function egloss_cont_new_with_real_v(key, val) result(cont)
    implicit none
    real, dimension(:), intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont=egloss_cont_new_with_egloss(key,list_new(.item. val))
  end function egloss_cont_new_with_real_v

  function egloss_cont_new_with_double_v(key, val) result(cont)
    implicit none
    double precision, dimension(:), intent(in) :: val
    character(len = *), intent(in) :: key
    type(eglossary_container) :: cont
    cont=egloss_cont_new_with_egloss(key,list_new(.item. val))
  end function egloss_cont_new_with_double_v
  
  !implementation of the .item. interface
  elemental function item_char(val) result(elem)
    implicit none
    character(len=*), intent(in) :: val
    type(list_container) :: elem
    elem%val(1:max_field_length)=val
  end function item_char

  elemental function item_integer(val) result(elem)
    implicit none
    integer, intent(in) :: val
    type(list_container) :: elem
    elem%val(1:max_field_length)=yaml_toa(val)
  end function item_integer

  elemental function item_real(val) result(elem)
    implicit none
    real, intent(in) :: val
    type(list_container) :: elem
    elem%val(1:max_field_length)=yaml_toa(val)
  end function item_real

  elemental function item_double(val) result(elem)
    implicit none
    double precision, intent(in) :: val
    type(list_container) :: elem
    elem%val(1:max_field_length)=yaml_toa(val)
  end function item_double

  function item_egloss(val) result(elem)
    implicit none
    type(eglossary), pointer, intent(in) :: val
    type(list_container) :: elem
    elem%egloss=>val
  end function item_egloss

  !implementation of the list_new interface
  function list_new(eglosses)
    implicit none
    type(list_container), dimension(:) :: eglosses
    type(eglossary), pointer :: list_new
    !local variables
    integer :: i_st,n_st
    type(eglossary), pointer :: egloss_tmp

    !initialize eglossionary
    call egloss_init(egloss_tmp)
    !TODO: need to implement add methods
    n_st=size(eglosses)
    do i_st=1,n_st
       if (associated(eglosses(i_st)%egloss)) then
          call add(egloss_tmp,eglosses(i_st)%egloss)
       else if (len_trim(eglosses(i_st)%val) > 0) then
          call add(egloss_tmp,eglosses(i_st)%val)
       end if
    end do

    list_new => egloss_tmp
  end function list_new
  !Create a list from several optional values (string or egloss).
  function list_new_elems(egloss0, egloss1, egloss2, egloss3, egloss4, egloss5, egloss6, egloss7, egloss8, egloss9)
    implicit none
    type(list_container), intent(in) :: egloss0
    type(list_container), intent(in), optional :: egloss1, egloss2, egloss3, egloss4
    type(list_container), intent(in), optional :: egloss5, egloss6, egloss7, egloss8, egloss9
    type(eglossary), pointer :: list_new_elems
    !local variables
    type(eglossary), pointer :: egloss_tmp

    !initialize eglossionary
    call egloss_init(egloss_tmp)
    call fill(egloss_tmp, egloss0)
    if (present(egloss1)) call fill(egloss_tmp, egloss1)
    if (present(egloss2)) call fill(egloss_tmp, egloss2)
    if (present(egloss3)) call fill(egloss_tmp, egloss3)
    if (present(egloss4)) call fill(egloss_tmp, egloss4)
    if (present(egloss5)) call fill(egloss_tmp, egloss5)
    if (present(egloss6)) call fill(egloss_tmp, egloss6)
    if (present(egloss7)) call fill(egloss_tmp, egloss7)
    if (present(egloss8)) call fill(egloss_tmp, egloss8)
    if (present(egloss9)) call fill(egloss_tmp, egloss9)
    list_new_elems => egloss_tmp
  contains
    subroutine fill(egloss, elem)
      implicit none
      type(list_container), intent(in) :: elem
      type(eglossary), pointer :: egloss

      if (associated(elem%egloss)) then
         call add(egloss, elem%egloss)
      else if (len_trim(elem%val) > 0) then
         call add(egloss, elem%val)
      end if
    end subroutine fill
  end function list_new_elems

  !implementation of the add iterface
  !> Add to a list
  subroutine add_char(egloss,val, last_item_ptr)
    implicit none
    type(eglossary), pointer :: egloss
    character(len=*), intent(in) :: val
    type(eglossary), pointer, optional :: last_item_ptr
    include 'simple_egloss_add-incFl.f90'
  end subroutine add_char
  subroutine add_egloss(egloss,val, last_item_ptr)
     implicit none
     type(eglossary), pointer :: egloss
     type(eglossary), pointer :: val
     type(eglossary), pointer, optional :: last_item_ptr
     include 'simple_egloss_add-incFl.f90'
   end subroutine add_egloss
   subroutine add_integer(egloss,val, last_item_ptr)
     implicit none
     type(eglossary), pointer :: egloss
     integer(kind=4), intent(in) :: val
     type(eglossary), pointer, optional :: last_item_ptr
     include 'simple_egloss_add-incFl.f90'
   end subroutine add_integer
   subroutine add_real(egloss,val, last_item_ptr)
     implicit none
     type(eglossary), pointer :: egloss
     real, intent(in) :: val
     type(eglossary), pointer, optional :: last_item_ptr
     include 'simple_egloss_add-incFl.f90'
   end subroutine add_real
   subroutine add_double(egloss,val, last_item_ptr)
     implicit none
     type(eglossary), pointer :: egloss
     double precision, intent(in) :: val
     type(eglossary), pointer, optional :: last_item_ptr
     include 'simple_egloss_add-incFl.f90'
   end subroutine add_double
   subroutine add_long(egloss,val, last_item_ptr)
     implicit none
     type(eglossary), pointer :: egloss
     integer(kind=8), intent(in) :: val
     type(eglossary), pointer, optional :: last_item_ptr
     include 'simple_egloss_add-incFl.f90'
   end subroutine add_long
   subroutine add_log(egloss,val, last_item_ptr)
     implicit none
     type(eglossary), pointer :: egloss
     logical, intent(in) :: val
     type(eglossary), pointer, optional :: last_item_ptr
     include 'simple_egloss_add-incFl.f90'
   end subroutine add_log

  !iterator
  function egloss_iter(egloss)
    implicit none
    type(eglossary),pointer,intent(in) :: egloss
    type(eglossary),pointer :: egloss_iter
    if (associated(egloss)) then
       egloss_iter=>egloss%child
    else
       nullify(egloss_iter)
    end if
  end function egloss_iter

  function egloss_iter_lc(list)
    implicit none
    type(list_container),intent(in) :: list
    type(eglossary), pointer :: egloss_iter_lc
    egloss_iter_lc => egloss_iter(list%egloss)
  end function egloss_iter_lc

  function egloss_next(egloss)
    implicit none
    type(eglossary),pointer,intent(in) :: egloss
    type(eglossary),pointer :: egloss_next
    !TODO: need to implement the method
    if (associated(egloss)) then
       if (associated(egloss%parent)) then
          egloss_next=>egloss%next
       else
          egloss_next=>egloss%child
       end if
    else
       nullify(egloss_next)
    end if
  end function egloss_next

  !getters
  !character
  recursive subroutine get_value(val,egloss)
    implicit none
    character(len=*), intent(out) :: val
    type(eglossary), intent(in) :: egloss
    !TODO: implement the checker code for the method
    val(1:len(val))=' '
    call file_strcpy(src=egloss%data%value,dest=val)
  end subroutine get_value
  !integer
  recursive subroutine get_integer(ival,egloss)
    use simple_yaml_strings, only: is_atoi
    implicit none
    integer(kind=4), intent(out) :: ival
    type(eglossary), intent(in) :: egloss
    !TODO: need to implement the method
  end subroutine get_integer
  !real
  subroutine get_real(rval,egloss)
    implicit none
    real(kind=4), intent(out) :: rval
    type(eglossary), intent(in) :: egloss
    !local variablles
    integer :: err
    double precision :: dval
    character(len=max_field_length) :: val
    val=egloss
    call read_fraction_string(val,dval,err)
    rval=real(dval)
    return
  end subroutine get_real
  !double
  subroutine get_double(dval,egloss)
    implicit none
    real(kind=8), intent(out) :: dval
    type(eglossary), intent(in) :: egloss
    !local variables
    integer :: ierror
    character(len=max_field_length) :: val
     
    !take value
    val=egloss
    !look at conversion
    !TODO: implement the read_fraction string to handle the decimal
    !call read_fraction_string(val, dval, ierror)

    if (file_err(ierror/=0,'Value '//val,err_id=EGLOSS_CONVERSION_ERROR)) return

    return
  end subroutine get_double

  
  !include the error_handling common interfaced methods and convoluted methods
  include 'simple_error_handling-incFl.f90'
   
end module simple_eglossary
