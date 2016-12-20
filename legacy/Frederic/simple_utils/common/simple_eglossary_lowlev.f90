! ============================================================================
! Name        : simple_glossary_lowlev
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 29th of June 2016
! Description : Module to define a glossary of low level fucntionalities,
!             : called by the simple_eglossary module. This includes the
!             : initialisers for the glossary.
! ============================================================================
!
module simple_eglossary_lowlev

  implicit none
  integer, parameter, public :: max_field_length=256    !Max field length
  character(len = max_field_length), parameter :: NOT_A_VALUE="__not_a_value__"
  character(len=max_field_length), parameter :: TYPE_EGLOSS='__egloss__'
  character(len=max_field_length), parameter :: TYPE_LIST='__list__'
  integer, private :: negloss=0
  integer, private :: negloss_max=0
  integer, private :: nfolders=0       !Number of libraries
  integer, private :: nfolders_max=0   !Max number of libraries allocated
  integer, parameter, private :: nfolder_size=10000 !Folder size of pre-allocated eglossary
  !storage type
  type, public :: storage
     sequence
     integer :: item
     integer :: nitems
     integer :: nelems
     character(len=max_field_length) :: key
     character(len=max_field_length) :: value
  end type storage
  !eglossary type 
  type, public :: eglossary
     type(storage) :: data
     type(eglossary), pointer :: parent => null()
     type(eglossary), pointer :: next => null()
     type(eglossary), pointer :: child => null()
     type(eglossary), pointer :: previous => null()
  end type eglossary
  
  !Database book is a workspace of pre-allocated eglossaries to manage
  !eglossary creation
  type,private :: database_book
     integer,dimension(:),pointer :: free_folders=>null()
     integer :: nfree
     type(eglossary), dimension(:),pointer :: files=>null()
     type(database_book),pointer :: next=>null()
     type(database_book),pointer :: previous=>null()
  end type database_book
  
  !Database type
  type(database_book), pointer, private :: library=>null()
  
  !Operator to access and create a key in the eglossary
  interface operator(//)
     module procedure get_child_ptr, get_list_ptr
  end interface
  !private methods
  private :: allocate_library,allocate_file,deallocate_library, deallocate_file
  !public declaration of the routine and functions
  public :: operator(//)
  public :: egloss_init, egloss_free
  public :: define_parent
  !getters
  public :: get_egloss_from_the_key
  !Handlers
  public :: egloss_len,egloss_size
  
contains

  !routine to return the value of the glossary egloss
  pure function egloss_value(egloss)
    implicit none
    type(eglossary), pointer, intent(in) :: egloss
    character(len=max_field_length) :: egloss_value

    if (associated(egloss)) then 
       !call check_key(egloss)
       if (associated(egloss%child)) then
          if (egloss%data%nitems > 0) then
             egloss_value=TYPE_LIST
          else if (egloss%data%nelems > 0) then
             egloss_value=TYPE_EGLOSS
          else
             egloss_value= NOT_A_VALUE !illegal condition
          end if
       else if (trim(egloss%data%value) == NOT_A_VALUE) then
          egloss_value=repeat(' ',len(egloss_value))
       else
          egloss_value=egloss%data%value
       end if
    else
       egloss_value=repeat(' ',len(egloss_value))
    end if


  end function egloss_value
    
  ! Return the size of the eglossionary
  pure function egloss_size(egloss)
    implicit none
    type(eglossary), intent(in), pointer :: egloss
    integer :: egloss_size
    
    if (associated(egloss)) then
       egloss_size=egloss%data%nelems
    else
       egloss_size=-1
    end if
  end function egloss_size
  !function that returns the v
  pure function egloss_item(egloss)
    type(eglossary),pointer, intent(in) :: egloss
    integer :: egloss_item

    if (associated(egloss)) then 
       egloss_item=egloss%data%item
    else
       egloss_item=-1
    end if

  end function egloss_item
  !Returns key value of eglossary
  pure function egloss_key(egloss)
    type(eglossary), pointer, intent(in) :: egloss
    character(len=max_field_length) :: egloss_key
    
    if (associated(egloss)) then 
       egloss_key=egloss%data%key
    else
       egloss_key=repeat(' ',len(egloss_key))
    end if
  end function egloss_key
  
  !this fucntion returns the key
  pure function name_is(egloss,name)
    implicit none
    type(eglossary),pointer, intent(in) :: egloss
    character(len=*), intent(in) :: name
    logical :: name_is
    name_is=.false.    
    if (trim(name) == trim(egloss%data%key)) then
       name_is=.true.
    else if (associated(egloss%child)) then
       !most likely egloss is a item list
       name_is=(trim(name) == trim(egloss%child%data%key))
       !write(*,*)'here',name_is,trim(egloss%child%data%value),'ag',trim(name)
    else
       name_is=(trim(name) == trim(egloss%data%value))
    end if
  end function name_is

  !allocate file
  function allocate_file() result(egloss)
    implicit none
    type(eglossary),pointer :: egloss
    !local variables
    integer :: ifolder
    type(database_book), pointer :: lib_tmp

    !create the place for a folder if none
    if (.not. associated(library%free_folders)) then
       call allocate_library(library)
       nfolders=nfolders+1
       nfolders_max=max(nfolders_max,nfolders)
    end if
    ! We use the saved free place.
    ifolder = library%free_folders(nfolder_size - library%nfree + 1)
    call egloss_nullify(library%files(ifolder))
    egloss=>library%files(ifolder)
    !and immediately put it unavailable
    library%nfree = library%nfree - 1
    library%free_folders(nfolder_size - library%nfree)=0

    do while (library%nfree == 0)
       !then the library which is used is the following one
       if (.not. associated(library%next)) then
          allocate(library%next) !the pointer
          lib_tmp => library%next
          lib_tmp%previous => library
          call allocate_library(library%next)
          nfolders=nfolders+1
          nfolders_max=max(nfolders_max,nfolders)
       end if
       library => library%next
    end do
  end function allocate_file
    
  !allocating the library and data initialising the database_book data structure
  pure subroutine allocate_library(library)
    implicit none
    type(database_book), intent(inout):: library
    !local variables
    integer :: ifolder
    allocate(library%free_folders(nfolder_size))
    allocate(library%files(nfolder_size))
    nullify(library%next)
    library%nfree = nfolder_size
    do ifolder =1,nfolder_size
       library%free_folders(ifolder) = ifolder
    end do
    return
  end subroutine allocate_library

  !Cleaners
  !deallocate library
  recursive subroutine deallocate_library
    implicit none
    type(eglossary),pointer :: egloss
    !local variables
    !TODO: need to implement the deallocation method
    return
  end subroutine deallocate_library
  !deallocate file
  subroutine deallocate_file(file)
    implicit none
    integer(kind=8) :: file
    !local variables
    !TODO: need to implement the deallocation method
  contains
  end subroutine deallocate_file

  !initialisers
  
  !eglossary initiliaser
  subroutine egloss_init(egloss)
    implicit none
    type(eglossary), pointer :: egloss
    if (nfolder_size == 0 ) then
       allocate(egloss)
       call egloss_nullify(egloss)
    else
       !allocate the files here
       if (.not. associated(library)) allocate(library)
       egloss=>allocate_file()
    end if
    negloss = negloss + 1
    negloss_max = max(negloss_max,negloss)
    return
  end subroutine egloss_init

  ! Retrieve the pointer to the item of the list.
  subroutine item_ptr_find(egloss,item,item_ptr)
    implicit none
    type(eglossary), intent(in), pointer :: egloss 
    integer, intent(in) :: item
    type(eglossary), pointer :: item_ptr

    item_ptr=>egloss
    find_item: do 
       if (item_ptr%data%item == item) exit find_item
       if (associated(item_ptr%next)) then
          item_ptr=>item_ptr%next
          cycle find_item
       end if
       if (.not. no_key(item_ptr)) then
          call init_next(item_ptr)
       end if
       call set_item(item_ptr,item)
       exit find_item
    end do find_item
  end subroutine item_ptr_find

  subroutine init_next(egloss)
    implicit none
    type(eglossary), pointer :: egloss

    call egloss_init(egloss%next)
    call define_brother(egloss,egloss%next) !chain the list in both directions
    if (associated(egloss%parent)) &
         call define_parent(egloss%parent,egloss%next)
    egloss=>egloss%next          
  end subroutine init_next

  
  !replaced with proper freer
  !  subroutine egloss_free()
  !    implicit none
  !    !TODO: insert the glossary initialisator, throw to yaml if required
  !    write(*,*) "need to free the glossary"
  !
  !    return
  !  end subroutine egloss_free

  !nullfier
  pure function storage_null() result(st)
    type(storage) :: st
    st%key=repeat(' ',max_field_length)
    st%value(1:max_field_length)=NOT_A_VALUE
    st%item = -1
    st%nitems=0
    st%nelems=0
  end function storage_null

  !nullifyer for the egloss pointer
  subroutine egloss_nullify(egloss)
    implicit none
    type(eglossary) :: egloss
    egloss%data = storage_null()
    nullify(egloss%parent,egloss%next,egloss%child,egloss%previous) 
    return
  end subroutine egloss_nullify

  !Return the length of the list
  pure function egloss_len(egloss)
    implicit none
    type(eglossary), intent(in), pointer :: egloss
    integer :: egloss_len
    
    if (associated(egloss)) then
       egloss_len=egloss%data%nitems
    else
       egloss_len=-1
    end if
  end function egloss_len
  !recursive routine to define the parent from the structure
  recursive subroutine define_parent(egloss,child)
    implicit none
    type(eglossary), target :: egloss
    type(eglossary) :: child
    !local variables
    type(eglossary),pointer :: iter
    child%parent=>egloss
    if (associated(child%next)) call define_parent(egloss,child%next)
!    return
  end subroutine define_parent
  !recursive routine to define the brother from the structure
  recursive subroutine define_brother(brother,egloss)
    implicit none
    type(eglossary), target :: brother
    type(eglossary) :: egloss
    egloss%previous => brother
!    return
  end subroutine define_brother
  
  !sets an elem in the tree 
  subroutine set_elem(egloss,key)
    use simple_yaml_strings
    implicit none
    type(eglossary), pointer :: egloss !!must be verified
    character(len=*), intent(in) :: key
    call file_strcpy(src=trim(key),dest=egloss%data%key)
    if (associated(egloss%parent)) then
       egloss%parent%data%nelems=egloss%parent%data%nelems+1
    else
       egloss%data%nelems=egloss%data%nelems+1
    end if
    write(*,*) 'set_elem out ',trim(key),egloss%data%nelems,egloss%parent%data%nelems
    return
  end subroutine set_elem

  !sets item in the tree 
  subroutine set_item(egloss,item)
    implicit none
    type(eglossary), pointer :: egloss !!must be verified
    integer, intent(in) :: item
    egloss%data%item=item
    if (associated(egloss%parent)) then
       if (len_trim(egloss%parent%data%value)>0)&
            egloss%parent%data%value=repeat(' ',max_field_length)
       egloss%parent%data%nitems=egloss%parent%data%nitems+1
    end if

    return
  end subroutine set_item
  
  !fucntion to get egloss from the key
  function get_egloss_from_the_key(egloss,key,create) result (egloss_ptr)
    implicit none
    !the root of the eglossary to start the search from
    type(eglossary), intent(in), pointer :: egloss
    character(len=*), intent(in) :: key
    logical, intent(in), optional :: create
    type(eglossary), pointer :: egloss_ptr
    !local variables
    logical :: crt
    type(eglossary), pointer :: iter

    crt = .false.
    if (present(create)) crt= create
    !iterator to find the key
    nullify(egloss_ptr)
    iter => egloss
    seek: do
       if (iter%data%key == trim(key) ) then
          egloss_ptr => iter
          exit seek
       else if (associated(iter%next)) then
          iter => iter%next
          cycle seek
       else
          exit seek
       end if
    end do seek
    if (crt) then
       if (no_key(iter)) then
          call set_elem(iter,key)
          egloss_ptr => iter
       end if
       !if not found create it
       if (.not. associated(egloss_ptr)) then
          call egloss_init(iter%next)
          call define_brother(iter,iter%next)!chains the list in both directions
          if (associated(iter%parent)) call define_parent(iter%parent,iter%next)
          call set_elem(iter%next,key)
          egloss_ptr => iter%next
       end if
    end if
    
  end function get_egloss_from_the_key

  !no_key function tests if the key is present
  pure function no_key(egloss)
    implicit none
    type(eglossary), intent(in) :: egloss
    logical :: no_key
    no_key = (len_trim(egloss%data%key) == 0 .and. egloss%data%item == -1) .and. associated(egloss%parent)
  end function no_key
    
  !getting the child pointer
  function get_child_ptr(egloss,key) result(sube_ptr)
    implicit none
    type(eglossary), intent(in), pointer :: egloss
    character(len=*),intent(in) :: key
    type(eglossary), pointer :: sube_ptr

    if (.not. associated(egloss) ) then
       nullify(sube_ptr)
       return
    end if

    if (associated(egloss%child)) then
       !TODO: implement the getter for egloss from the child
       sube_ptr => get_egloss_from_the_key(egloss%child,key,create=.true.)
    else
       call egloss_init(egloss%child)
       call define_parent(egloss,egloss%child)
       call set_elem(egloss%child,key)
       sube_ptr=>egloss%child
    end if
    
  end function get_child_ptr
  !gettin the list pointer
  function get_list_ptr(egloss,item) result(sube_ptr)
    type(eglossary), intent(in), pointer :: egloss
    integer, intent(in) :: item
    type(eglossary), pointer :: sube_ptr

    if (.not. associated(egloss) ) then
       nullify(sube_ptr)
       return
    end if
    call clean_subegloss(egloss)
    if (associated(egloss%child)) then
       !TODO: check if implementation of item_ptr_find is correct
       call item_ptr_find(egloss%child,item,sube_ptr)
    else
       call egloss_init(egloss%child)
       call define_parent(egloss,egloss%child)
       call set_item(egloss%child,item)
       sube_ptr=>egloss%child
    end if
    
  end function get_list_ptr

  !cleaner for sub egloss
  subroutine clean_subegloss(egloss)
    implicit none
    type(eglossary), pointer :: egloss
    if (associated(egloss)) then
       if (egloss%data%nelems > 0) then
          call egloss_free(egloss%child)
          egloss%data%nelems = 0
       end if
    end if
    return
  end subroutine clean_subegloss

  !freeing the eglossary
  subroutine egloss_free(egloss)
    implicit none
    type(eglossary), pointer :: egloss
    if (associated(egloss)) then
       call egloss_free_(egloss)
       nullify(egloss)
    end if
  contains
    !recursive sunbroutine to free the tree
    recursive subroutine egloss_free_(egloss)
      implicit none
      type(eglossary), pointer :: egloss
      !local variables
      type(eglossary),pointer :: egloss_iter,child,current
      egloss_iter => egloss
      do while (associated(egloss_iter))
         child => egloss_iter%child
         current => egloss_iter
         egloss_iter => egloss_iter%next
         !destroy current
         if(nfolder_size == 0) then
            deallocate(current)
         else
            call egloss_destroy(egloss)
         end if
         nullify(current)
         !destroy the children
         if (associated(child)) then
            call egloss_free_(child)
         end if
      end do
    end subroutine egloss_free_

  end subroutine egloss_free

  !subroutine to destroy the eglossary
  subroutine egloss_destroy(egloss)
    implicit none
    type(eglossary), intent(in) :: egloss
    !local variables

    !TODO: need to implement the destroy method

  end subroutine egloss_destroy
  
end module simple_eglossary_lowlev
